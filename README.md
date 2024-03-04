# Testing Random Polynomials for Polynomial Conservation Laws

This algorithm was developed as the output of my internship during the summer of 2023, supervised by Dr. Rahkooy, University of Oxford. 

Written in Python, depending primarily on SymPy, this algorithm finds the syzygy module bases of an inputted ideal (which comes from ODEbase), generates a random linear combination of the syzygies, checks if the curl of the new syzygy module bases is equal to zero, if it is equal to zero then this is syzygy is integrated. The final integration is then called the conservation law (sometimes known as first integrals).

To use this algorithm, one needs to define your model as seen in the experiment code and then run the genpolconslaw function with the variables defined as one can observe in the experiment code.  

The next step in this project is to write a bash script which converts the models from ODEbase to the type that is used by the algorithm.


ODEbase: https://odebase.org/table/?num_species_range=&num_parameters_range=&num_constraints_range=&deficiency_range=&is_rational=&is_polynomial=True&is_mass_action=
