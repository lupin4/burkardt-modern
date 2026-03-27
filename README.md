# burkardt-modern

### John Burkardt's Computational Geometry — Modernized to Fortran 2018

53 libraries upgraded from Fortran 90 reference code to production-grade Fortran 2018 with parallelism, vectorization, and modern type safety.

---

## What Changed

| Before (Burkardt original) | After (modernized) |
|---|---|
| `integer ( kind = 4 )` | `integer(int32)` via `iso_fortran_env` |
| `real ( kind = 8 )` | `real(real64)` via `iso_fortran_env` |
| Serial `do` loops | `do concurrent` where parallelizable |
| No purity annotations | `pure` / `elemental` on all eligible routines |
| No `intent` on many params | Explicit `intent(in/out/inout)` everywhere |
| `implicit none` per routine | `implicit none` via module-level |
| No `contiguous` | `contiguous` on assumed-shape arrays |
| No OpenMP | `!$omp parallel do` on hot loops |
| `kind = 4/8` magic numbers | Named constants from `iso_fortran_env` |
| No `bind(C)` | `bind(C, name="...")` on all public routines |
| Standalone subroutines | Organized into modules |
| No error handling convention | `info` integer status codes (0 = success) |

---

## Libraries

### Core Geometry
| Library | Lines | Description |
|---------|-------|-------------|
| **geometry** | 41,100 | 500+ routines for 2D/3D/ND computational geometry |
| **dutch** | 6,180 | Computational geometry (Dutch school algorithms) |

### Triangulation & Meshing
| Library | Lines | Description |
|---------|-------|-------------|
| **geompack3** | 44,900 | Barry Joe — 2D/3D/ND Delaunay + Voronoi |
| **geompack2** | 14,016 | Extended Delaunay with hole handling |
| **geompack** | 4,042 | Original geompack |
| **triangulation** | 14,586 | 2D triangulation, order 3 and 6 |
| **triangulation_quality** | 2,633 | Mesh quality metrics |
| **polygon_triangulate** | 2,282 | Polygon ear-clipping triangulation |
| **delaunay_lmap_2d** | 3,996 | Delaunay with linear map |
| **table_delaunay** | 3,173 | Delaunay from point sets |

### Sphere Operations
| Library | Lines | Description |
|---------|-------|-------------|
| **stripack** | 8,526 | Robert Renka — spherical Delaunay/Voronoi |
| **sphere_grid** | 3,316 | Multiple sphere grid strategies |
| **sphere_quad** | 3,538 | Quadrature on sphere surfaces |
| **sphere_fibonacci_grid** | 504 | Fibonacci spiral sphere grid |
| **sphere_cubed_grid** | 1,190 | Cube-projected sphere grid |
| **sphere_llt_grid** | 770 | Lat/lon triangle grid |
| **sphere_llq_grid** | 736 | Lat/lon quad grid |
| **sphere_delaunay** | 1,587 | Spherical Delaunay triangulation |
| **sphere_voronoi** | 1,826 | Spherical Voronoi diagram |
| **sphere_stereograph** | 1,452 | Stereographic projection |

### Grid Generation
| Library | Lines | Description |
|---------|-------|-------------|
| **ball_grid** | 554 | Grid points inside 3D ball |
| **cube_grid** | 496 | Grid points inside cube |
| **disk_grid** | 558 | Grid points inside 2D disk |
| **ellipse_grid** | 587 | Grid points inside ellipse |
| **ellipsoid_grid** | 626 | Grid points inside 3D ellipsoid |
| **hex_grid_angle** | 1,430 | Hexagonal grid generation |
| **hypercube_grid** | 549 | Grid inside M-dimensional hypercube |
| **line_grid** | 209 | Grid on 1D line segment |
| **polygon_grid** | 597 | Grid inside polygon |
| **pyramid_grid** | 637 | Grid inside pyramid |
| **simplex_grid** | 1,116 | Grid inside M-dimensional simplex |
| **tetrahedron_grid** | 158 | Grid inside tetrahedron |
| **triangle_grid** | 158 | Grid inside triangle |
| **wedge_grid** | 458 | Grid inside wedge |
| **circle_arc_grid** | 367 | Grid along circular arc |

### Tetrahedral Mesh
| Library | Lines | Description |
|---------|-------|-------------|
| **tet_mesh** | 6,233 | Tet mesh operations |
| **tet_mesh_boundary** | 3,086 | Boundary extraction |
| **tet_mesh_quality** | 3,761 | Quality metrics |
| **tetrahedron_grid** | 569 | Grid inside tetrahedron |
| **tetrahedron_properties** | 2,971 | Circumsphere, insphere, volume, angles |

### Polygon & Polyhedron
| Library | Lines | Description |
|---------|-------|-------------|
| **polygon_properties** | 3,973 | Area, centroid, containment, diameter |
| **polygon_integrals** | 431 | Exact polygon integrals |
| **polygon_grid** | 597 | Grid inside polygon |
| **circle_segment** | 3,020 | Circle segment computations |

### Quadrilateral Mesh
| Library | Lines | Description |
|---------|-------|-------------|
| **quad_mesh** | 5,198 | Quad mesh operations |
| **quad_mesh_rcm** | 4,597 | RCM bandwidth reduction |

### Other
| Library | Lines | Description |
|---------|-------|-------------|
| **felippa** | 6,456 | Quadrature on every geometric primitive |
| **hypersphere_properties** | 1,864 | N-dimensional hypersphere ops |
| **simplex_coordinates** | 557 | Regular simplex vertices in M-D |
| **cvt** | 3,030 | Centroidal Voronoi Tessellation |
| **bezier_surface** | 3,285 | Bezier surface evaluation |
| **tri_surface_io** | 2,192 | Triangle surface I/O |
| **naca** | 541 | NACA airfoil geometry |

---

## Modernization Status

| Library | Original | Modern | Status |
|---------|----------|--------|--------|
| geometry | `src/original/` | `src/modern/` | Pending |
| ... | ... | ... | ... |

Each library progresses: **Original** -> **Typed** -> **Modularized** -> **Parallelized** -> **Complete**

---

## Build

```bash
# Build all modernized libraries
make all

# Build specific library
make geometry

# Run tests (original vs modern output comparison)
make test

# Benchmark original vs modern
make bench
```

---

## Directory Structure

```
burkardt-modern/
+-- src/
|   +-- original/          Unmodified Burkardt sources (reference)
|   +-- modern/            Modernized Fortran 2018 sources
|   +-- tests/             Regression tests (original == modern)
+-- benchmarks/            Performance comparisons
+-- scripts/               Modernization tooling
+-- Makefile
+-- README.md
```

---

## Modernization Checklist (per library)

- [ ] Replace `kind = 4/8` with `iso_fortran_env` named constants
- [ ] Add `implicit none` at module level
- [ ] Mark `pure` / `elemental` where eligible
- [ ] Add explicit `intent(in/out/inout)` on all parameters
- [ ] Add `contiguous` on assumed-shape array arguments
- [ ] Convert eligible `do` loops to `do concurrent`
- [ ] Add `bind(C, name="...")` on public routines
- [ ] Wrap in modules with `private` default + explicit `public`
- [ ] Add OpenMP directives on parallelizable hot loops
- [ ] Fix column-major access patterns where row-major is used
- [ ] Add `info` status code return where appropriate
- [ ] Regression test: modern output matches original output

---

## Why

The Burkardt collection is the gold standard for computational geometry in Fortran — correct, documented, battle-tested across decades. But it's written in Fortran 90 style with no parallelism, no purity annotations, and magic `kind` numbers.

This repo modernizes the code while preserving correctness. Every modernized routine is regression-tested against the original. The result: the same trusted algorithms, running faster, with compile-time safety guarantees.

**Upstream for:** [forGeo](https://github.com/forKernels/forGeo) (forKernels computational geometry module)

---

## Credits

Original code by **John Burkardt** (University of South Carolina)
https://people.math.sc.edu/Burkardt/f_src/

Original license: **GNU LGPL**
Modernization: Same license.
