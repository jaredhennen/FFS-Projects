These files represent a few Python and IDL programs I wrote during the course of my research at the University of Minnesota, as well as a Jupyter Notebook for a problem found on "The Riddler" at fivethirtyeight.com.

Riddler 9_29.ipynb - Solution to a problem on the Riddler. A baseball league consists of 3 teams, one which walks for 40% of at bats (Moonwalkers), one which hits doubles for 20% of at bats (Doublooners), and one which hits home runs for 10% of at bats (Taters), all teams strike out otherwise. Find the probability of each team winning a season. Here I used a monte carlo simulation approach to simulate 10,000 seasons and determine the probabilities, while also visualizing the distribution of wins using the seaborn library.

RunZ-Scan.py - Python program with a GUI interface used to output a waveform while collecting photon counts for FFS experiments. The user then selects points of interest to perform longer data collection, which the program then collects automatically. This program greatly enhanced the accuracy of data collection while simultaneously reducing the time for each experiment.

createbootstrapmsqarray.pro - IDL program used to perform bootstrapped tsMSQ analysis as described in Hennen et al 2019, Analytical Biochemistry.

hoppingrandomwalk.pro - IDL program for a single iteration of a random walk experiment which was then repeated. This was used to test the feasability of a proposed experiment involving the scanning of the microscope's observation volume. Goal was to differentiate between simple diffusion and "hopping" diffusion which involves stochastic immobilization (but with an identical effective diffusion coefficient).

jhmsq_define.pro - Object oriented IDL program used for some basic MSQ analysis. Typically used for on the fly analysis during experiments but is used as a precurser to bootstrap analysis.
