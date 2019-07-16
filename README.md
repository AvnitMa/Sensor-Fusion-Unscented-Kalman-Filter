# **Sensor Fusion: Unscented Kalman Filter**

In this project I implemented an unscented Kalman filter using the CTRV motion model. I am using the same bicycle simulation data set from the extended Kalman filter project. 

## Reflection

### What is Unscented Kalman Filter?

* All Kalman filters have the same three steps:

    1. Initialization
    2. Prediction
    3. Update

* A standard Kalman filter can only handle linear equations. 
Both the extended Kalman filter and the unscented Kalman filter allow you to use non-linear equations.
The difference between EKF and UKF is how they handle non-linear equations.
The extended Kalman filter can give particularly poor performance.

* For example: in a car, if you asume that the velocity and direction are constant, you can use the extended Kalman filter. 
However, since the road as many turns, the model will predict turning vehicles incorrectly, becuase the position estimation whould tend to result outside the actually driven circle.

![](http://i.imgur.com/X5II4S0.png)

* The unscented Kalman filter uses a sampling technique known as the unscented transform to pick a sigma points around the mean. 

* These sigma points are then propagated through the non-linear functions to get a new mean and covariance estimatations.

* This filter can more accurately estimates the true mean and covariance.

* An additional benefit is that we don't need to calculate Jacobians (which is a complicated function).


#### Ouput video:

[![Unc](http://i.imgur.com/SXFYVOJ.png)](https://vimeo.com/229897312 "Unc")
