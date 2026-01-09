//! Spatial indexing for efficient neighbor search.
//!
//! This module provides a simple grid-based spatial index for finding
//! atoms within a cutoff radius.

use std::collections::HashMap;

/// Grid-based spatial index for 3D point queries.
///
/// Divides space into uniform cubic cells and stores atom indices
/// in each cell. Supports efficient range queries within a cutoff radius.
#[derive(Debug)]
pub struct SpatialGrid {
    /// Inverse cell size for fast coordinate-to-cell conversion.
    inv_cell_size: f64,
    /// Map from cell coordinates to atom indices.
    cells: HashMap<(i32, i32, i32), Vec<usize>>,
}

impl SpatialGrid {
    /// Creates a new spatial grid with the given cell size.
    ///
    /// # Arguments
    ///
    /// * `cell_size` — Size of each cubic cell (typically the cutoff radius)
    ///
    /// # Panics
    ///
    /// Panics if `cell_size <= 0.0`.
    pub fn new(cell_size: f64) -> Self {
        assert!(cell_size > 0.0, "Cell size must be positive");
        Self {
            inv_cell_size: 1.0 / cell_size,
            cells: HashMap::new(),
        }
    }

    /// Creates a spatial grid and populates it with atom positions.
    ///
    /// # Arguments
    ///
    /// * `positions` — Slice of 3D positions [x, y, z] in Ångströms
    /// * `cell_size` — Size of each cubic cell (typically the cutoff radius)
    pub fn from_positions(positions: &[[f64; 3]], cell_size: f64) -> Self {
        let mut grid = Self::new(cell_size);
        for (idx, pos) in positions.iter().enumerate() {
            grid.insert(idx, *pos);
        }
        grid
    }

    /// Computes the cell coordinates for a given position.
    fn cell_coords(&self, pos: [f64; 3]) -> (i32, i32, i32) {
        (
            (pos[0] * self.inv_cell_size).floor() as i32,
            (pos[1] * self.inv_cell_size).floor() as i32,
            (pos[2] * self.inv_cell_size).floor() as i32,
        )
    }

    /// Inserts an atom index at the given position.
    ///
    /// # Arguments
    ///
    /// * `idx` — Atom index to store
    /// * `pos` — 3D position of the atom
    pub fn insert(&mut self, idx: usize, pos: [f64; 3]) {
        let cell = self.cell_coords(pos);
        self.cells.entry(cell).or_default().push(idx);
    }

    /// Finds all atom indices within the cutoff radius of any query point.
    ///
    /// This is optimized for the case where we have multiple query points
    /// (e.g., all atoms of a ligand) and want to find all environment atoms
    /// within range of any query point.
    ///
    /// # Arguments
    ///
    /// * `queries` — Slice of query positions
    /// * `positions` — Full position array for distance calculations
    /// * `cutoff` — Maximum distance to include
    ///
    /// # Returns
    ///
    /// Sorted, deduplicated vector of atom indices within cutoff of any query.
    pub fn query_radius_multi(
        &self,
        queries: &[[f64; 3]],
        positions: &[[f64; 3]],
        cutoff: f64,
    ) -> Vec<usize> {
        let cutoff_sq = cutoff * cutoff;

        let mut cells_to_search = std::collections::HashSet::new();
        for query in queries {
            let (cx, cy, cz) = self.cell_coords(*query);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        cells_to_search.insert((cx + dx, cy + dy, cz + dz));
                    }
                }
            }
        }

        let mut result_set = std::collections::HashSet::new();
        for cell in cells_to_search {
            if let Some(indices) = self.cells.get(&cell) {
                for &idx in indices {
                    let pos = positions[idx];
                    for query in queries {
                        let dist_sq = (pos[0] - query[0]).powi(2)
                            + (pos[1] - query[1]).powi(2)
                            + (pos[2] - query[2]).powi(2);
                        if dist_sq <= cutoff_sq {
                            result_set.insert(idx);
                            break;
                        }
                    }
                }
            }
        }

        let mut results: Vec<_> = result_set.into_iter().collect();
        results.sort_unstable();
        results
    }

    /// Finds all atom indices within the cutoff radius of a query point.
    ///
    /// This is a single-point query variant, primarily used in tests.
    /// For batch queries, use [`query_radius_multi`](Self::query_radius_multi).
    ///
    /// # Arguments
    ///
    /// * `query` — Query position [x, y, z]
    /// * `positions` — Full position array for distance calculations
    /// * `cutoff` — Maximum distance to include
    ///
    /// # Returns
    ///
    /// Vector of atom indices within the cutoff radius.
    #[cfg(test)]
    pub fn query_radius(&self, query: [f64; 3], positions: &[[f64; 3]], cutoff: f64) -> Vec<usize> {
        let cutoff_sq = cutoff * cutoff;
        let (cx, cy, cz) = self.cell_coords(query);

        let mut results = Vec::new();

        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let cell = (cx + dx, cy + dy, cz + dz);
                    if let Some(indices) = self.cells.get(&cell) {
                        for &idx in indices {
                            let pos = positions[idx];
                            let dist_sq = (pos[0] - query[0]).powi(2)
                                + (pos[1] - query[1]).powi(2)
                                + (pos[2] - query[2]).powi(2);
                            if dist_sq <= cutoff_sq {
                                results.push(idx);
                            }
                        }
                    }
                }
            }
        }

        results
    }
}
