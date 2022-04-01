---
date: '2022-04-01T13:55:18.499Z'
docname: Examples/developer_examples/dev_example_topology_creation
images: {}
path: /examples-developer-examples-dev-example-topology-creation
title: Topology creation and modifications with PyGromosTools
---

# Topology creation and modifications with PyGromosTools

PyGromosTools offers a wide variaty of tools and functions to assist in the creation and modification of specialized topologies.

In the following notebook we will demonstrate some of the functions and theire usage. These functions and workflows should work with all supported forcfields, but we will use OpenForceField due to the fact that it is independant of make_top and a Gromos++ installation.

## Automatic System generation and combining of topologies

In the first part of this notebbok we demonstrate some automatic tools to create and modifie topologies.

After the automatic creation of this topology for a Cyclohexan molecule we could try to make a liquid topology. Therefore we want to have multiple copies of our topology in one single topology.

With the “+=” operator we can add multiple topologies to a single topology. Or multiply the topology in a for loop. This alows costum topologies for liquids. After the loop we can reasign the variable to our Gromos_System and check if the additional molecules were added.

We can also add different molecules in our system. For example we could add a Toluene to our Cyclohexane system.

If many topologies should be added to a single topology the PyGromosTools version of com_top can be used which directly implements multipliers. Instead of for loops.

## Manual topology modifications

In this second part of the notebook we demonstrate some functions to manually tweek topologies.

Like with any othe PyGromosTools file class we have the option to brute force rewrite blocks. However, this often requires advanced knowlede of the block structure. Except for a few simple cases like the TITLE block.

For more complexe blocks PyGromosTools provides many functions to modify topologies. For example we could simply add a new soluteatom to the system. This would be relevant if you want to manually add a ion with very specific Nonbonded values to you system. These tools are vary powrful, but don’t do any sanity checks on the values or compatibility to your forcefields.

<script type="application/vnd.jupyter.widget-state+json">
{"state": {"59752c5e4fd246eab4c71b2d31954721": {"model_name": "LayoutModel", "model_module": "@jupyter-widgets/base", "model_module_version": "1.2.0", "state": {"_model_module": "@jupyter-widgets/base", "_model_module_version": "1.2.0", "_model_name": "LayoutModel", "_view_count": null, "_view_module": "@jupyter-widgets/base", "_view_module_version": "1.2.0", "_view_name": "LayoutView", "align_content": null, "align_items": null, "align_self": null, "border": null, "bottom": null, "display": null, "flex": null, "flex_flow": null, "grid_area": null, "grid_auto_columns": null, "grid_auto_flow": null, "grid_auto_rows": null, "grid_column": null, "grid_gap": null, "grid_row": null, "grid_template_areas": null, "grid_template_columns": null, "grid_template_rows": null, "height": null, "justify_content": null, "justify_items": null, "left": null, "margin": null, "max_height": null, "max_width": null, "min_height": null, "min_width": null, "object_fit": null, "object_position": null, "order": null, "overflow": null, "overflow_x": null, "overflow_y": null, "padding": null, "right": null, "top": null, "visibility": null, "width": null}}, "d2538164463747dfb4394761fd9fdaa6": {"model_name": "ColormakerRegistryModel", "model_module": "nglview-js-widgets", "model_module_version": "3.0.1", "state": {"_dom_classes": [], "_model_module": "nglview-js-widgets", "_model_module_version": "3.0.1", "_model_name": "ColormakerRegistryModel", "_msg_ar": [], "_msg_q": [], "_ready": false, "_view_count": null, "_view_module": "nglview-js-widgets", "_view_module_version": "3.0.1", "_view_name": "ColormakerRegistryView", "layout": "IPY_MODEL_59752c5e4fd246eab4c71b2d31954721"}}}, "version_major": 2, "version_minor": 0}
</script>
