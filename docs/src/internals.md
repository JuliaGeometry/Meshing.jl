# Internals

The internals of Meshing have been optimized for performance and to be generic.
This is some brief documentation on the basic internal functions.

## Common

```@docs
Meshing._DEFAULT_SAMPLES
Meshing._get_cubeindex
```

## Marching Cubes

```@docs
Meshing._mc_create_triangles!
Meshing._mc_unique_triangles!
Meshing.find_vertices_interp
Meshing.vertex_interp
Meshing.mc_vert_points
```


## Marching Tetrahedra

```
Voxel corner and edge indexing conventions

        Z
        |

        5------5------6              Extra edges not drawn
       /|            /|              -----------
      8 |           6 |              - face diagonals
     /  9          /  10                - 13: 1 to 3
    8------7------7   |                 - 14: 1 to 8
    |   |         |   |                 - 15: 1 to 6
    |   1------1--|---2  -- Y           - 16: 5 to 7
    12 /          11 /                  - 17: 2 to 7
    | 4           | 2                   - 18: 4 to 7
    |/            |/                 - body diagonal
    4------3------3                     - 19: 1 to 7

  /
 X
```

```@docs
Meshing._correct_vertices!
Meshing.procVox
Meshing.voxEdgeId
Meshing.getVertId
Meshing.vertPos
Meshing.vertId
Meshing.tetIx
```

## Surface Nets