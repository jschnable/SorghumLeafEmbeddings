LeafWebScore — Leaf Disease Severity Scoring Tool
==================================================

A simple tool for scoring leaf disease severity from images.
Runs locally in your web browser. No installation required beyond Python 3.


REQUIREMENTS
------------
- Python 3 (pre-installed on Mac and Linux; install from python.org on Windows)
- A web browser (Chrome, Firefox, Safari, Edge, etc.)


SETUP
-----
1. Place your image folders inside the "ImageProjects" directory.
   Each folder becomes a project in the scoring tool.

   Example:
     ImageProjects/
       AAMU/
         image1.jpg
         image2.jpg
         image3.jpg
       FVSU/
         photo_a.jpg
         photo_b.jpg

   - Folder names become the project names shown in the dropdown menu.
   - Any common image format works (JPG, PNG, GIF, BMP, TIFF, WebP).
   - Images can be nested in subdirectories if needed — the tool will find them.
   - The tool ignores hidden files and macOS system files (__MACOSX, .DS_Store, etc.).

2. Start the server:

     python3 server.py

   This will automatically open your web browser to the scoring page.
   If it doesn't, open http://localhost:8000 manually.


SCORING
-------
1. Enter your name and click "Start Scoring".
2. Select a project from the dropdown at the top.
3. For each image, click a score button (1.0 to 7.0) then click Submit.
   - Or use keyboard shortcuts for speed:
       1-7     Select score (e.g. press 4 for 4.0)
       .       Toggle the .5 (e.g. 4.0 becomes 4.5)
       Enter   Submit the score
       S       Skip the image (for blurry/unreadable images)
       B       Go back to the previous image
4. Optionally type a comment before submitting.
5. Your progress is saved automatically. You can close and come back later —
   just enter the same name to resume where you left off.

Image navigation:
  - Scroll wheel to zoom in/out
  - Click and drag to pan
  - Double-click to reset zoom


OUTPUT
------
Scores are saved to "scores.csv" in this directory. Each row is one score:

  project,image,username,score,timestamp,comment

- Scores are saved immediately (crash-safe).
- If you re-score an image, both entries are kept. The latest one is the valid score.
- "skip" in the score column means the scorer skipped that image.


TIPS
----
- Multiple people can score at the same time by opening the URL in different
  browsers (or on different computers on the same network).
- Each scorer sees images in a different random order to reduce bias.
- The same scorer always gets the same order, so progress is consistent.
- To add a new project, just add a new folder to ImageProjects/ and restart
  the server.
