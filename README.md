# FlowLB – How to Run

## HOW TO RUN (Linux Terminal or VS Code Terminal using Makefile – CPU)

1. The Makefiles are already provided for both CPU:

   * `Makefile`

2. You may need to adjust the relative folder locations of the input files in `main.cpp`:

   * `params.txt`
   * `settings.txt`

   (You can copy settings and params files OR create a pointer link to any example in examples/ by ln -s examples/Case1_parametersLB_2D.txt config/params.txt)

3. Commands for clean and build are as follows.(If you are running the program from 1_FlowLB_2D_3D_MPI/ you first need to change the address of settings.txt in main from "../../config/settings.txt" to "config/settings.txt". The same for params.txt):

   * `make -f Makefile clean`
   * `make -f Makefile`

4. Run the program:

   ```
   mpirun -np 5 ./FlowLB_CPU
   ```

---

## HOW TO RUN (Code::Blocks – CPU Only)

1. Open `FlowLB_3D_MPI.cbp` in Code::Blocks.
2. Input files are located in:

   * `config/settings.txt`
   * `config/params.txt`
     Adjust them as needed.
3. Compile and build the code.
4. Run the executable from the `Debug/` or `Release/` folder.
   If you need 9 processors for flow, use 10 processors (one processor is reserved for particles):
   ```bash
      mpirun -n 10 FlowLB_2D_3D
   ```
5. Validate with the non-parallel version located in:

   ```
   This_version_is_validated_with/
   ```

---

## NOTES

* This version works only for:

  ```
  num_processors_for_particles_only = 1
  ```
* To run 2D cases:

  * Set `Z_Periodic`
  * Use `Nz = 1`
  * Do not use any thermal boundary conditions
  * Particle types are limited to **circle** or **ellipse**
* In the current version:

  * Only the **MASTER** process (`processor_id = 0`) handles particles
  * All other processors handle fluid
  * The MASTER receives fluid data for interpolation and distributes forces after spreading

