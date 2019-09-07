
"""
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
"""
const voxCrnrPos = (Point(0, 0, 0),
                Point(0, 1, 0),
                Point(1, 1, 0),
                Point(1, 0, 0),
                Point(0, 0, 1),
                Point(0, 1, 1),
                Point(1, 1, 1),
                Point(1, 0, 1))
# the voxel IDs at either end of the tetrahedra edges, by edge ID
const voxEdgeCrnrs = ((1, 2),
                (2, 3),
                (4, 3),
                (1, 4),
                (5, 6),
                (6, 7),
                (8, 7),
                (5, 8),
                (1, 5),
                (2, 6),
                (3, 7),
                (4, 8),
                (1, 3),
                (1, 8),
                (1, 6),
                (5, 7),
                (2, 7),
                (4, 7),
                (1, 7))

# direction codes:
# 0 => +x, 1 => +y, 2 => +z,
# 3 => +xy, 4 => +xz, 5 => +yz, 6 => +xyz
const voxEdgeDir = (1,0,1,0,1,0,1,0,2,2,2,2,3,4,5,3,4,5,6)

# For a pair of corner IDs, the edge ID joining them
# 0 denotes a pair with no edge
const voxEdgeIx = ((0,1,13,4,9,15,19,14),
                (1,0,2,0,0,10,17,0),
                (13,2,0,3,0,0,11,0),
                (4,0,3,0,0,0,18,12),
                (9,0,0,0,0,5,16,8),
                (15,10,0,0,5,0,6,0),
                (19,17,11,18,16,6,0,7),
                (14,0,0,12,8,0,7,0))

# voxel corners that comprise each of the six tetrahedra
const subTets = ((1,3,2,7),
                (1,8,4,7),
                (1,4,3,7),
                (1,2,6,7),
                (1,5,8,7),
                (1,6,5,7))
# tetrahedron corners for each edge (indices 1-4)
const tetEdgeCrnrs = ((1,2),
                (2,3),
                (1,3),
                (1,4),
                (2,4),
                (3,4))

# triangle cases for a given tetrahedron edge code
const tetTri = ((1,3,4,0,0,0),
            (1,5,2,0,0,0),
            (3,5,2,3,4,5),
            (2,6,3,0,0,0),
            (1,6,4,1,2,6),
            (1,5,6,1,6,3),
            (4,5,6,0,0,0),
            (4,6,5,0,0,0),
            (1,6,5,1,3,6),
            (1,4,6,1,6,2),
            (2,3,6,0,0,0),
            (3,2,5,3,5,4),
            (1,2,5,0,0,0),
            (1,4,3,0,0,0))
