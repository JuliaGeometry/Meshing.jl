
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
const voxCrnrPos = ((0, 0, 0),
                    (0, 1, 0),
                    (1, 1, 0),
                    (1, 0, 0),
                    (0, 0, 1),
                    (0, 1, 1),
                    (1, 1, 1),
                    (1, 0, 1))

# the voxel IDs at either end of the tetrahedra edges, by edge ID
const voxEdgeCrnrs = ((0x01, 0x02),
                      (0x02, 0x03),
                      (0x04, 0x03),
                      (0x01, 0x04),
                      (0x05, 0x06),
                      (0x06, 0x07),
                      (0x08, 0x07),
                      (0x05, 0x08),
                      (0x01, 0x05),
                      (0x02, 0x06),
                      (0x03, 0x07),
                      (0x04, 0x08),
                      (0x01, 0x03),
                      (0x01, 0x08),
                      (0x01, 0x06),
                      (0x05, 0x07),
                      (0x02, 0x07),
                      (0x04, 0x07),
                      (0x01, 0x07))

# direction codes:
# 0 => +x, 1 => +y, 2 => +z,
# 3 => +xy, 4 => +xz, 5 => +yz, 6 => +xyz
const voxEdgeDir = (0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0x02, 0x02, 0x02, 0x02, 0x03, 0x04, 0x05, 0x03, 0x04, 0x05, 0x06)

# For a pair of corner IDs, the edge ID joining them
# 0 denotes a pair with no edge
const voxEdgeIx = ((0x00, 0x01, 0x0d, 0x04, 0x09, 0x0f, 0x13, 0x0e),
                   (0x01, 0x00, 0x02, 0x00, 0x00, 0x0a, 0x11, 0x00),
                   (0x0d, 0x02, 0x00, 0x03, 0x00, 0x00, 0x0b, 0x00),
                   (0x04, 0x00, 0x03, 0x00, 0x00, 0x00, 0x12, 0x0c),
                   (0x09, 0x00, 0x00, 0x00, 0x00, 0x05, 0x10, 0x08),
                   (0x0f, 0x0a, 0x00, 0x00, 0x05, 0x00, 0x06, 0x00),
                   (0x13, 0x11, 0x0b, 0x12, 0x10, 0x06, 0x00, 0x07),
                   (0x0e, 0x00, 0x00, 0x0c, 0x08, 0x00, 0x07, 0x00))

# voxel corners that comprise each of the six tetrahedra
const subTets = ((0x01, 0x03, 0x02, 0x07),
                 (0x01, 0x08, 0x04, 0x07),
                 (0x01, 0x04, 0x03, 0x07),
                 (0x01, 0x02, 0x06, 0x07),
                 (0x01, 0x05, 0x08, 0x07),
                 (0x01, 0x06, 0x05, 0x07))

# voxel corners that comprise each of the six tetrahedra, masking cubeindex
const subTetsMask = ((0x04, 0x02),
                     (0x80, 0x08),
                     (0x08, 0x04),
                     (0x02, 0x20),
                     (0x10, 0x80),
                     (0x20, 0x10))

# tetrahedron corners for each edge (indices 1-4)
const tetEdgeCrnrs = ((0x01, 0x02),
                      (0x02, 0x03),
                      (0x01, 0x03),
                      (0x01, 0x04),
                      (0x02, 0x04),
                      (0x03, 0x04))

# triangle cases for a given tetrahedron edge code
const tetTri = ((0x01, 0x03, 0x04, 0x00, 0x00, 0x00),
                (0x01, 0x05, 0x02, 0x00, 0x00, 0x00),
                (0x03, 0x05, 0x02, 0x03, 0x04, 0x05),
                (0x02, 0x06, 0x03, 0x00, 0x00, 0x00),
                (0x01, 0x06, 0x04, 0x01, 0x02, 0x06),
                (0x01, 0x05, 0x06, 0x01, 0x06, 0x03),
                (0x04, 0x05, 0x06, 0x00, 0x00, 0x00),
                (0x04, 0x06, 0x05, 0x00, 0x00, 0x00),
                (0x01, 0x06, 0x05, 0x01, 0x03, 0x06),
                (0x01, 0x04, 0x06, 0x01, 0x06, 0x02),
                (0x02, 0x03, 0x06, 0x00, 0x00, 0x00),
                (0x03, 0x02, 0x05, 0x03, 0x05, 0x04),
                (0x01, 0x02, 0x05, 0x00, 0x00, 0x00),
                (0x01, 0x04, 0x03, 0x00, 0x00, 0x00))
