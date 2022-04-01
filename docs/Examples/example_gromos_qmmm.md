---
date: '2022-04-01T13:55:18.499Z'
docname: Examples/example_gromos_qmmm
images: {}
path: /examples-example-gromos-qmmm
title: QM/MM in Gromos
---

# QM/MM in Gromos

## Support for QMMM functionality in `GROMOS` input files

This notebook demonstrates support of `PyGromosTools` for QM/MM functionality.

[https://github.com/rinikerlab/PyGromosTools/blob/qmmm/examples/example_gromos_qmmm.ipynb](https://github.com/rinikerlab/PyGromosTools/blob/qmmm/examples/example_gromos_qmmm.ipynb) (part of the `qmmm` branch and soon to be merged to `release3`)

Author: Felix Pultar

Features include:


* QM/MM blocks in `imd` files


* QM/MM specification files


* Running QM/MM simulations

### Load an `imd` file containing a QMMM block

Simple demonstration of how to handle `.imd` files.

### Print out different sections of the QMMM block

Print out selected parameters from the `QMMM` block or also the `TITLE` block.

### Change a block value and print again

Just change values of the `QMMM` block like with other `PyGromosTools` blocks.

## Directly manipulate a QMMM specification file

The QMMM object allows to directly interact with QM/MM specification files. Future releases of `PyGromosTools` will also support generation of `QMMM` files from coordinate files (`.cnf`, `.xyz`, `.pdb`).

## Print out and change some blocks in the QMMM specification file

### Title block

The `QMMM` specification file can be handled like any other `GROMOS` file.

### QMZONE block

Print out the `QMZONE` section that defines which atoms will be treated quantum-mechanically.

### QMUNIT block

Print out the `QMUNIT` block that defines some unit conversions between the MD engine and the QM software.

### XTBELEMENTS block

Print and update the `XTBELEMENTS` block

### A helper function that returns all QM engines specified in the QM/MM specification file

There is also a sanity check in the constructor of `QMMM` to see if you did not accidentally add more than one QM engine.

### Store your QMMM specification file with all your other simulation files in a `Gromos_System` object

## Run QM/MM Simulations

QM/MM calculations are possible with a special in-house build of `GROMOS`.

<script type="application/vnd.jupyter.widget-state+json">
{"state": {"817a4de7320249a5b64bc982f12f62a6": {"model_name": "LayoutModel", "model_module": "@jupyter-widgets/base", "model_module_version": "1.2.0", "state": {"_model_module": "@jupyter-widgets/base", "_model_module_version": "1.2.0", "_model_name": "LayoutModel", "_view_count": null, "_view_module": "@jupyter-widgets/base", "_view_module_version": "1.2.0", "_view_name": "LayoutView", "align_content": null, "align_items": null, "align_self": null, "border": null, "bottom": null, "display": null, "flex": null, "flex_flow": null, "grid_area": null, "grid_auto_columns": null, "grid_auto_flow": null, "grid_auto_rows": null, "grid_column": null, "grid_gap": null, "grid_row": null, "grid_template_areas": null, "grid_template_columns": null, "grid_template_rows": null, "height": null, "justify_content": null, "justify_items": null, "left": null, "margin": null, "max_height": null, "max_width": null, "min_height": null, "min_width": null, "object_fit": null, "object_position": null, "order": null, "overflow": null, "overflow_x": null, "overflow_y": null, "padding": null, "right": null, "top": null, "visibility": null, "width": null}}, "0fac11cf325d489d8dc7a6fc87300ad1": {"model_name": "ColormakerRegistryModel", "model_module": "nglview-js-widgets", "model_module_version": "3.0.1", "state": {"_dom_classes": [], "_model_module": "nglview-js-widgets", "_model_module_version": "3.0.1", "_model_name": "ColormakerRegistryModel", "_msg_ar": [], "_msg_q": [], "_ready": false, "_view_count": null, "_view_module": "nglview-js-widgets", "_view_module_version": "3.0.1", "_view_name": "ColormakerRegistryView", "layout": "IPY_MODEL_817a4de7320249a5b64bc982f12f62a6"}}}, "version_major": 2, "version_minor": 0}
</script>
