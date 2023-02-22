# Autopilot-Design-with-PSO

In this project we are trying to find best time domain performance from the autopilot where we obey the frequency
domain requirements. Obviously, this is an optimization problem and solved by "particle swarm optimization(PSO)".
The algorithm that we used in PSO solver is modified by us such that; it restart itselves every 1/5 th of the 
total iteration number. By doing that we are trying to reduce the "probability of local minima convergence problem"
and results show that we are correct.

There is a primitive GUI inside the code, which can be improved and used more effectively. Robustness criteria and time 
domain performance requirements are taken from a paper which is the basement of this project. I will give the link below so 
you can go and read the entire paper about the topic.

https://doi.org/10.1016/j.ifacol.2016.03.161
Siddhardha Kedarisetty
