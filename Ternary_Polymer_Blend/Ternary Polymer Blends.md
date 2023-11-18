# Ternary Polymer Blends



## Postprocessing

Use the following postprocessing steps to generate the desired figure:

1. Apply the `ExtractBlock` filter to $a$ and $b$ fields separately, extracting only the variable and not the mesh. 
2. For each field, select either the `Reds`, `Blues`, or `Greens` color map. 
3. Select the `Invert the transfer function` option in the color map editor. 
4. Select `Enable Opacity Mapping For Surfaces` which makes the color map transparent when $c_{i} = 0$
5. 

