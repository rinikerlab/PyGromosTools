---
date: '2022-04-01T13:55:18.499Z'
docname: Examples/example_gromos_on_HPC_Clusters
images: {}
path: /examples-example-gromos-on-hpc-clusters
title: Usage of the PyGromosTools on HPC-cluster
---

# Usage of the PyGromosTools on HPC-cluster

Note that you need support for the new `PyGromosTools` submission system shipped with `relase3`.

## Initialization

Set up path to binaries and instantiate an initial `Gromos_System` object

## Energy minimization

Take advantage of the new submission system: convenience functions such as `emin` take in a parametrized `Gromos_System` object, automatically set up the calculation and return a new `Gromos_System` object with updated file paths.

### Print out the new file paths

Note how the new `Gromos_System` object has new file paths associated with it.

### Visualize the minimized system

**Note:** The new `Gromos_System` object has an .imd file associated (from the energy minimization). Using the object for a subsequent calculation requires the .imd to be reset (as shown below) or a new `Imd` object to be passed to the function `md`. Otherwise `md` will be based on the previous (emin) input file.

## Equilibration followed by subsequent production runs

Again, note how the new `Gromos_System` object has new file paths associated with it.

### Visualize the last configuration

<script type="application/vnd.jupyter.widget-state+json">
{"state": {"1773c0fb6da14520996b886c56b1c191": {"model_name": "LayoutModel", "model_module": "@jupyter-widgets/base", "model_module_version": "1.2.0", "state": {"_model_module": "@jupyter-widgets/base", "_model_module_version": "1.2.0", "_model_name": "LayoutModel", "_view_count": null, "_view_module": "@jupyter-widgets/base", "_view_module_version": "1.2.0", "_view_name": "LayoutView", "align_content": null, "align_items": null, "align_self": null, "border": null, "bottom": null, "display": null, "flex": null, "flex_flow": null, "grid_area": null, "grid_auto_columns": null, "grid_auto_flow": null, "grid_auto_rows": null, "grid_column": null, "grid_gap": null, "grid_row": null, "grid_template_areas": null, "grid_template_columns": null, "grid_template_rows": null, "height": null, "justify_content": null, "justify_items": null, "left": null, "margin": null, "max_height": null, "max_width": null, "min_height": null, "min_width": null, "object_fit": null, "object_position": null, "order": null, "overflow": null, "overflow_x": null, "overflow_y": null, "padding": null, "right": null, "top": null, "visibility": null, "width": null}}, "85ba40b172b04866827c9e153d63b315": {"model_name": "ColormakerRegistryModel", "model_module": "nglview-js-widgets", "model_module_version": "3.0.1", "state": {"_dom_classes": [], "_model_module": "nglview-js-widgets", "_model_module_version": "3.0.1", "_model_name": "ColormakerRegistryModel", "_msg_ar": [], "_msg_q": [], "_ready": false, "_view_count": null, "_view_module": "nglview-js-widgets", "_view_module_version": "3.0.1", "_view_name": "ColormakerRegistryView", "layout": "IPY_MODEL_1773c0fb6da14520996b886c56b1c191"}}}, "version_major": 2, "version_minor": 0}
</script>
