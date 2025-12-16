# CBET
Code for the CBET paper: parallel and sequential versions in Rust and C++. Compile with `cargo build --release` in the Rust case and `make` in the C++ case.

The constants governing the size of the mesh, the number of rays used per laser beam, beam intensity, and other details are in `src/consts.rs` and `consts.hpp` respectively. This simplified simulation of the phenomenon is limited to two laser beams in a 2-dimensional mesh with a linear electron gradient.
