# PHY571

## Group 8:
### Cazenave Nicolas, Lam Kon Seng Léo, Brunel Clément

Numerical physics project on Bird flocking and Swarm intelligence modelisation

Here is the organisation of the different folders:
    
        - design_code: 
the first draft of the classes we want to build later on
not useful to understand the project, as we left it behind early on

        - our_library:
library composed of three classes:
    - Simple_model: the implementation of the code of Viscek's article
    - Complex_model: the implementation of a more complex model, which takes into accounts other interactions such as trying to avoid birds too close.
    - Predator_model: the implementation of a predator in the flock.

        - production_code: 
allows to do the animations for the different models
allows to plot the different figures of Viscke's article, both with his model and our more complex model

        - results:
here the data (that is generated and plotted thanks to production_code) is saved.
