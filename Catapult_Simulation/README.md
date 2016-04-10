# README  --  Catapult Simulation

This directory contains Matlab code to simulate a simple torsion catapult. This tutorial is designed to show event detection on a two-phase simulation problem, with some new plotting and animation tricks as well.

### Simulation:

__Part One:__ The catapult accelerates the projectile, powered by a torsional spring. The simulation terminates when the catapult arm reaches a prescribed angle.

__Part Two:__ The projectile is now released and flies through the air, under the influence of quadratic drag and gravity. The simulation terminates the when the projectile reaches the ground.

### Model:

__Catapult:__ The arm of the catapult is a slendar rod, and the projectile (sitting at the tip of the arm) is a point-mass. The arm is accelerated using a torsion spring.

__Projectile:__ Point mass, projectile motion with quadratic drag.

__Ground:__ The ground is modeled as a sum of sine and cosine functions, making the final event detection problem more interesting. 

