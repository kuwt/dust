**Purpose:**

This repo is for the coupling with Dymore. There is no difference between this code and DUST, except logging of lifting line data.  Therefore, cloning from DUST directly should provide the same or better functionality. The only thing added by this repo is the test examples that is for the coupling with Dymore. They are locate at

./examples/dymore_wing

**Tutorial:**
1. install preCICE, DUST and Dymore
2. Make sure a DUST case is corresponded to a Dymore case.  
3. Go to the directory of a DUST case. Run dust_pre. Then run dust. You should be able to see DUST is waiting for preCICE connection
4. Go to the directory of a Dymore case. Run Dymore with the (.dym). You should be able to see Dymore is waiting for preCICE connection.
5. The order of steps 2,3 can be interchanged. 
6. In util/tool_for_PlottingAndModelGeneration, there are tools for plotting the logging added for this coupling.

**Correspondence betwen DUST case and Dymore case**
DUST                                  <->       Dymore
1. elliptic_wing or rectangular wing  <->       bo105/wing
2. rotor_2b                           <->       rotor_2b
3. rotor_4b                           <->       rotor_4b or bo105adv_test_dust
