#References


* **Volonteri, M. 2015 (1509.09027)**

    Brief review of issue related to BHs in galaxy mergers.

* **Springel, V. 2004 (0411108)**

    Basics of our black hole growth model.

* **Springel, V. 2002 (0206393)**
    
    Star formation mode.

* **Sijacki, D. 2007 (0705.2238)**

    Extension of the model with radio mode.

* **Sijacki, D. 2011 (1008.3313)**

    A first study on the effects of gravitational recoil.

* **Petts, J. A. 2015 (1509.07871)**

    A recent idea to devise a dynamical friction formula for use in N-body models.

* **Tremmel, M. 2015 (1501.07609)**

    Another paper in this direction, which could be a starting point for our own numerical work.

* **Kahn, F. 2016 (1604.00015)**

    Swift coalsecence of binary SMBH in cosmological mergers. A zoom in simulation is performed
    using the GASOLINE code. Here a couple of merging galaxies are taken from a Cosmological 
    simulation of the project ARGO. The ratio of the background particles (the most massive) and 
    the BHs is approximately 1:100, so dynamical friction is well accounted for by the numerical 
    integration. Once the binary is hardened and DF becomes inefficient, the central region of 
    the remnant is further integrated using N-body direct sum, where three bodies processes are 
    dominant. Finally, at the scale of a few parsecs and below, post newtonian terms are included 
    in order to account for gravitational wave emissions. The conclusions are that the inclusion 
    of PN terms shortens the decay time significantly compared with previous works where mergers 
    of isolated galaxies have been studied. **As a future reference, it would be interesting to 
    reproduce the results of these paper using our semianalytical approach.**
    
* **Tamburello, V. 2016 (1603.00021)**

    They study the influence of clumps of star forming gas on the decay times of binary SMBHs at 
    redshifts z~1-3. This happends either by fragmentation of the disk due to instabilities or by 
    minor merger events. Contrary to what is commonly believed, the presence of gas rich 
    environments can delay the decay time of a companion BH rather than speeding it up, as giant
    molecular clouds can strongly influence the orbit of the BH. They use high resolution 
    simulations, where the mass ratio of the background particles and the BH is roughly 1:1000, so
    dynamical friction is well accounted for by integration. They use the Chandreskhar formula, but
    for analytical estimatives of decays times, not for integration.
    
* **Blecha, L. 2015 (1508.01524)**

    They study effects of recoiling black holes on spin aligments. They take mergers from Illustris
    simulation and resimulate (integrate) the trajectory of the BHs using the computed potential 
    and an extra Chandrasekhar term that accounts for the dynamical friction. DF is introduced this
    way because the BH integration during the simulation is not self-consistent, and a wake in the
    background medium is not produced, which in the end, is what causes the DF effect.