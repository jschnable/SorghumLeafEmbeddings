#!/usr/bin/env python3
"""
LeafWebScore — Zero-dependency Python server for leaf disease severity scoring.

Usage:
    python3 server.py

Place image directories inside ImageProjects/:
    ImageProjects/
        AAMU/
            image1.jpg
            image2.jpg
        FVSU/
            image1.jpg
        UNL/
            image1.jpg

Each subdirectory becomes a project in the dropdown.
Serves index.html and persists scores to scores.csv.
"""

import csv
import hashlib
import json
import mimetypes
import os
import random
import socket
import sys
import threading
import webbrowser
from datetime import datetime
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path
from urllib.parse import unquote

# ── Configuration ──────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
PORT = 8000
SCORE_FILE = SCRIPT_DIR / "scores.csv"
IMAGES_DIR = SCRIPT_DIR / "ImageProjects"
IMAGE_EXTENSIONS = {".jpg", ".jpeg", ".png", ".gif", ".bmp", ".tiff", ".tif", ".webp"}
CSV_LOCK = threading.Lock()

# ── Image discovery ────────────────────────────────────────────────────────────

def discover_projects():
    """Scan ImageProjects/ for subdirectories containing images.

    Robust to common setup mistakes:
    - Recursively finds images even if nested in subdirectories
    - Skips hidden files/dirs (starting with . or __MACOSX)
    - Accepts any common image extension regardless of case
    - Handles spaces and special characters in names
    - Skips empty directories silently
    """
    projects = {}

    if not IMAGES_DIR.exists():
        return projects

    for entry in sorted(IMAGES_DIR.iterdir()):
        if not entry.is_dir():
            continue
        # Skip hidden/system directories
        if entry.name.startswith(".") or entry.name.startswith("__"):
            continue

        # Recursively find all image files in this project directory
        images = []
        for f in entry.rglob("*"):
            if not f.is_file():
                continue
            if f.name.startswith(".") or f.name.startswith("__"):
                continue
            if f.suffix.lower() in IMAGE_EXTENSIONS:
                # Store path relative to the project dir so we can serve it
                images.append(str(f.relative_to(entry)))

        if not images:
            print(f"  WARNING: '{entry.name}' has no image files — skipping")
            continue

        images.sort()
        project_name = entry.name
        projects[project_name] = {
            "path": entry,
            "images": images,
        }

    return projects


def get_scored_images(project, username):
    """Return set of image filenames scored by this user for this project."""
    scored = set()
    if not SCORE_FILE.exists():
        return scored
    with CSV_LOCK:
        with open(SCORE_FILE, "r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row["project"] == project and row["username"] == username:
                    scored.add(row["image"])
    return scored


def get_user_order(images, username):
    """Deterministic per-user shuffle of image list."""
    seed = int(hashlib.sha256(username.encode()).hexdigest(), 16)
    rng = random.Random(seed)
    shuffled = list(images)
    rng.shuffle(shuffled)
    return shuffled


def append_score(project, image, username, score, comment=""):
    """Append a single score row to CSV (crash-safe: flush immediately)."""
    with CSV_LOCK:
        file_exists = SCORE_FILE.exists()
        with open(SCORE_FILE, "a", newline="") as f:
            writer = csv.writer(f)
            if not file_exists:
                writer.writerow(["project", "image", "username", "score", "timestamp", "comment"])
            writer.writerow([project, image, username, score,
                             datetime.now().isoformat(timespec="seconds"), comment])
            f.flush()
            os.fsync(f.fileno())

# ── HTTP Handler ───────────────────────────────────────────────────────────────

PROJECTS = {}  # populated at startup


class Handler(SimpleHTTPRequestHandler):
    def log_message(self, format, *args):
        if args and isinstance(args[0], str) and args[0].startswith("GET /api"):
            return
        if args and isinstance(args[0], str) and args[0].startswith("POST"):
            return
        super().log_message(format, *args)

    def send_json(self, data, status=200):
        body = json.dumps(data).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", len(body))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        path = unquote(self.path)

        # Serve index.html at root
        if path == "/" or path == "/index.html":
            return self._serve_file(str(SCRIPT_DIR / "index.html"), "text/html")

        # API: list projects
        if path == "/api/projects":
            result = {
                "projects": {
                    name: len(p["images"]) for name, p in PROJECTS.items()
                }
            }
            return self.send_json(result)

        # API: progress
        if path.startswith("/api/progress/"):
            parts = path.split("/")
            if len(parts) >= 5:
                project = parts[3]
                username = "/".join(parts[4:])
                if project not in PROJECTS:
                    return self.send_json({"error": "Unknown project"}, 404)
                scored = get_scored_images(project, username)
                total = len(PROJECTS[project]["images"])
                return self.send_json({
                    "scored": len(scored),
                    "total": total,
                    "scored_images": sorted(scored)
                })

        # API: next unscored image
        if path.startswith("/api/next/"):
            parts = path.split("/")
            if len(parts) >= 5:
                project = parts[3]
                username = "/".join(parts[4:])
                if project not in PROJECTS:
                    return self.send_json({"error": "Unknown project"}, 404)
                scored = get_scored_images(project, username)
                order = get_user_order(PROJECTS[project]["images"], username)
                for img in order:
                    if img not in scored:
                        idx = order.index(img)
                        return self.send_json({
                            "image": img,
                            "index": idx + 1,
                            "total": len(order)
                        })
                return self.send_json({"image": None, "done": True})

        # API: get image at specific position (for go-back)
        if path.startswith("/api/image_at/"):
            parts = path.split("/")
            if len(parts) >= 6:
                project = parts[3]
                username = parts[4]
                try:
                    index = int(parts[5])
                except ValueError:
                    return self.send_json({"error": "Invalid index"}, 400)
                if project not in PROJECTS:
                    return self.send_json({"error": "Unknown project"}, 404)
                order = get_user_order(PROJECTS[project]["images"], username)
                if 1 <= index <= len(order):
                    return self.send_json({
                        "image": order[index - 1],
                        "index": index,
                        "total": len(order)
                    })
                return self.send_json({"error": "Index out of range"}, 400)

        # Serve images: /images/<project>/<relative_path>
        if path.startswith("/images/"):
            rest = path[len("/images/"):]
            # Split into project + remaining path
            slash = rest.find("/")
            if slash > 0:
                project = rest[:slash]
                rel_path = rest[slash + 1:]
                if project in PROJECTS:
                    img_path = (PROJECTS[project]["path"] / rel_path).resolve()
                    project_root = PROJECTS[project]["path"].resolve()
                    # Prevent path traversal outside project directory
                    if not str(img_path).startswith(str(project_root) + os.sep) and img_path != project_root:
                        self.send_error(403)
                        return
                    if img_path.exists() and img_path.is_file():
                        content_type = mimetypes.guess_type(str(img_path))[0] or "application/octet-stream"
                        return self._serve_file(str(img_path), content_type)
            self.send_error(404)
            return

        self.send_error(404)

    def do_POST(self):
        path = unquote(self.path)

        if path == "/api/score":
            length = int(self.headers.get("Content-Length", 0))
            body = json.loads(self.rfile.read(length))
            project = body.get("project")
            image = body.get("image")
            username = body.get("username")
            score = body.get("score")
            comment = body.get("comment", "")

            if not all([project, image, username, score is not None]):
                return self.send_json({"error": "Missing fields"}, 400)
            if project not in PROJECTS:
                return self.send_json({"error": "Unknown project"}, 404)

            append_score(project, image, username, score, comment)

            scored = get_scored_images(project, username)
            total = len(PROJECTS[project]["images"])
            return self.send_json({
                "ok": True,
                "scored": len(scored),
                "total": total
            })

        self.send_error(404)

    def _serve_file(self, filepath, content_type):
        try:
            with open(filepath, "rb") as f:
                data = f.read()
            self.send_response(200)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", len(data))
            self.end_headers()
            self.wfile.write(data)
        except FileNotFoundError:
            self.send_error(404)

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    global PROJECTS

    print("LeafWebScore — Leaf Disease Severity Scoring Tool")
    print("=" * 50)

    if not IMAGES_DIR.exists():
        IMAGES_DIR.mkdir()
        print(f"\nCreated '{IMAGES_DIR}/' directory.")
        print(f"Add project folders with images inside it, then restart.")
        print(f"Example:")
        print(f"  {IMAGES_DIR}/")
        print(f"    AAMU/")
        print(f"      image1.jpg")
        print(f"      image2.jpg")
        print(f"    FVSU/")
        print(f"      image1.jpg")
        sys.exit(0)

    print(f"\nScanning {IMAGES_DIR}/ for projects...")
    PROJECTS = discover_projects()

    if not PROJECTS:
        print(f"\nERROR: No projects with images found in '{IMAGES_DIR}/'.")
        print(f"Each subdirectory of '{IMAGES_DIR}/' should contain image files.")
        print(f"Currently found directories:")
        for entry in IMAGES_DIR.iterdir():
            if entry.is_dir() and not entry.name.startswith("."):
                print(f"  {entry.name}/  (no images found)")
        sys.exit(1)

    print(f"\nProjects loaded:")
    for name, info in PROJECTS.items():
        print(f"  {name}: {len(info['images'])} images")

    if SCORE_FILE.exists():
        with open(SCORE_FILE, "r") as f:
            line_count = sum(1 for _ in f) - 1  # minus header
        if line_count > 0:
            print(f"\nExisting scores: {line_count} entries in {SCORE_FILE}")

    # Try to find an available port
    server = None
    actual_port = PORT
    for port in range(PORT, PORT + 20):
        try:
            server = HTTPServer(("0.0.0.0", port), Handler)
            actual_port = port
            break
        except OSError:
            if port == PORT:
                print(f"\n  Port {port} is in use, trying next...")
            continue

    if server is None:
        print(f"\nERROR: Could not find an available port ({PORT}–{PORT + 19}).")
        print("Close other programs using these ports and try again.")
        sys.exit(1)

    local_url = f"http://localhost:{actual_port}"
    print(f"\nServer running at {local_url}")

    # Show LAN address for other users on the network
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("8.8.8.8", 80))
        lan_ip = s.getsockname()[0]
        s.close()
        print(f"Other computers on your network can connect to: http://{lan_ip}:{actual_port}")
    except Exception:
        pass

    print("Press Ctrl+C to stop.\n")

    try:
        webbrowser.open(local_url)
    except Exception:
        pass

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
        server.shutdown()


if __name__ == "__main__":
    main()
