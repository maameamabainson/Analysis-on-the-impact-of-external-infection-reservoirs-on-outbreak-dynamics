# Analysis-on-the-impact-of-external-infection-reservoirs-on-outbreak-dynamics

Diseases that complicate eradication efforts include those with external reservoirs ---  a source of the infectious agent other than the host population. This research focuses on diseases in this category and investigates how an external infection reservoir influences persistence and fadeout. We achieve this by estimating a key measure of disease impact, the expected time to outbreak fadeout, for two stochastic disease models that incorporate different representations of the reservoir. Our main methodology uses the Wentzell-Kramers-Brillouin (WKB) framework, which until now has primarily been applied to models without considering such external sources. We derive an asymptotic expression that explicitly accounts for the contribution of animal reservoirs to the expected fadeout time. For the model that includes environmental reservoir risk, we numerically estimate this quantity. Diseases with environmental reservoirs were found to produce more prolonged outbreaks than those with animal reservoirs. A policy effectiveness analysis suggests that with increasing risk of external infections, strategies favouring curative measures are better suited to diseases with environmental reservoirs, whereas those with greater allocation to preventive measures were more effective for diseases with animal reservoirs.

- animal_reservoir_model_verifying_analytical_result.m
  
This code verifies the analytical results for the mean return time and
the mean fadeout time with the eigenvalue approach to obtaining the same
quantities.

   Usage:
   
     R0 = animal_reservoir_model_verifying_analytical_result(kappa)
  
  Returns:
  
      R0:                     Basic reproduction number
      
  Required:
  
      kappa:                  Choice of kappa parameter, N*force of external infection.

  
- optimal_trajectory_computation_environmental_model.m

This code computes the optimal trajectory (solution of the equations of motion) for the SIS model with
environmental reservoir using either time truncation or time transformation.

Run with 'optimal_trajectory_plotting.m

Usage:

 [data, extinct_path, runtime] = optimal_trajectory_computation_environmental_model('CTs')
 
 [data, extinct_path, runtime] = optimal_trajectory_computation_environmental_model('CTr')


Returns:

 data:                  [varying parameter, R0, action, max_Hmt, endemic_distance, df_distance, path_length]      
 
Required:

  algorithm:             Algorithm key: 'CTs' for time transformation;
                                        'CTr' for time truncation


- optimal_trajectory_plotting.m

Underlying code used to plot optimal trajectories in Figure 6 in the paper.

Run with 'optimal_trajectory_computation_environmental_model.m'

Usage:

   tok = optimal_trajectory_plotting(extinct_path, end_eq, df_eq,T) 
   
Returns:

   tok:                    Computational time
      
Required:

 extinct_path:           Computed extinction optimal trajectory.
 
 end_eq:                 Corresponding endemic equilibrium (start point)
 
 df_eq:                  Corresponding disease-free equilibrium (end point)
 
 T:                      Absolute value of truncated/transformed time endpoint



                                        
