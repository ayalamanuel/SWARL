# Surface Wave-Aerodynamic Roughness Length (SWARL) Model for Air-Sea Interactions

The **SWARL model** estimates the aerodynamic roughness length $z_0$ for turbulent flow over moving ocean waves. SWARL is physics-based model and only requires the following inputs:  
1. Surface elevation map of the wavy surface.  
2. Airflow conditions, specifically the turbulent Reynolds number (or equivalent inputs like friction velocity and air viscosity).  

For details of the model derivation and validation, see **Ayala et al. (2024)**.

---

## Repository Contents  

### 1. **Surfaces**  
A folder containing 30 ocean wave surfaces:  
- 12 Monochromatic waves 
- 18 Multiscale waves generated using the JONSWAP or Pierson-Moskowitz spectrum models.  
Each surface is described in **Ayala et al. (2024)**.  

### 2. **Model Scripts**  
- **`zo_wave_model_.m`**:  
   MATLAB script to compute $z_0$ for any type of wavy surface.    
- **`solve_U.m`**:  
   A MATLAB function implementing the Newton-Raphson method to iteratively solve for the air velocity $U$, used to compute $\Lambda$.  

---

## Usage Instructions   
1. Clone or download this repository.
2. Open **`zo_wave_model.m`** script.
3. Select a wave surface type (monochromatic or multiscale)
4. Select a wave surface case from the available options:  
   - Monochromatic: W1, Z1, Z2, H1, H2, C1, C2, B1–B5.
   - Multiscale: Y1, Y2, D1–D4, J1–J4, R1–R4, S1–S4.
5. Run the **`zo_wave_model.m`** script.  
  
The script will compute the normalized surface roughness length scale.    

---

## Using Custom Ocean Wave Surfaces  

The SWARL model can analyze other wave surfaces as long as you provide:  
1. Surface elevation map (grid of elevation values).  
2. Friction velocity ($u_*$) and air viscosity ($\nu$) or the friction Reynolds number ($Re_\tau$).  

**Steps for Custom Surfaces:**  
1. Prepare 2 realizations of your wavy surface elevation data.  
2. Normalize the wave phase velocities ($c_x, c_y$) by the friction velocity ($u_*$).  
3. Use the provided scripts to calculate $z_0$.  

---

## Authors  

* **[Manuel Ayala](https://www.linkedin.com/in/manuelayalag/)** - *Johns Hopkins University*  
* **Dennice F. Gayme** - *Johns Hopkins University*  
* **Charles Meneveau** - *Johns Hopkins University*  

---

## Reference  

If you use this model, please cite:  
**Ayala, M., Gayme, D. F., & Meneveau, C. (2024). Surface Wave-Aerodynamic Roughness Length Model for Air-Sea Interactions. Submitted  to *Geophysical Research Letters*.**

---

