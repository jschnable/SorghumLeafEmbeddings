# Sorghum leaf imaging and embedding genetics

Code, data and supporting resources for **Embeddings of standardized sorghum leaf images capture genetic variation in disease response missed by human scores**.

## Supporting information

- [Collection protocol](supporting_info/collection_protocol/collection_protocol.pdf): illustrated instructions for leaf imaging and FieldBook data collection.
- [Imaging chamber](supporting_info/imaging_chamber/): printable chamber parts and a parts list.
- [LeafWebScore](supporting_info/LeafWebScore/ReadMe.txt): a local browser application for scoring leaf disease severity.

## Data and analyses

- [Analysis scripts and usage](scripts/README.md): analysis scripts and shared helpers in `scripts/`, with installation instructions, a single-image example, workflows and tests.
- [Provided data](data/provided/README.md): observations, field metadata, analysis settings and population-structure inputs.
- [External inputs](data/externalsourcerequired/README.md): access instructions for images, model weights, genotypes, expression data and annotations.
- [Generated datasets](data/generatable/README.md): output locations and their roles in the analyses.
- [Figures](figures/README.md): manuscript figures, plotting inputs and rendering instructions, with generation and assembly scripts in each figure's subdirectory under `figures/main/` or `figures/supplemental/`.
- [Shared figure inputs](data/figure_inputs/README.md): inputs used across candidate-gene panels.

To reproduce the analyses, start with the [script instructions](scripts/README.md) and [input inventory](data/provided/README.md).

## License

Unless otherwise noted, original content in this repository is licensed under [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/). See [LICENSE](LICENSE) for the full terms. Third-party materials retain their respective licenses.
