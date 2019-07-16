#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
    
  // first measurement
  x_ << 1, 1, 1, 1, 1;

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
    
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.8;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.7;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

 
  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;
    
  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

    
  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

    is_initialized_ = false;
   
    time_us_ =0;

    //set state dimension
    n_x_ = 5;
    
    //set augmented dimension
    n_aug_ = 7;
    
    n_z_ = 3;
    
    lambda_ = -2;
    
    weights_ = VectorXd(2*n_aug_+1);
    
    //create augmented mean vector
    x_aug_ = VectorXd(7);
    
    //create augmented state covariance
    P_aug_ = MatrixXd(7, 7);
    
    //create sigma point matrix
    Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    
    Zsig_ = MatrixXd::Zero(n_z_, 2 * n_aug_ + 1);

    R_laser_ = MatrixXd(2, 2);
    
    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
    0, 0.0225;
   
    P_ << 1, 0, 0, 0, 0,
    0, 1, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 1;
    
    H_laser_ = MatrixXd(2, 5);
    
    H_laser_ << 1, 0, 0, 0, 0,
    0, 1, 0, 0, 0;
    
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    
    double weight_0 = lambda_/(lambda_+n_aug_);
    weights_(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
        double weight = 0.5/(n_aug_+lambda_);
        weights_(i) = weight;
    }

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
 
    
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        
        
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            /**
             Convert radar from polar to cartesian coordinates and initialize state.
             */
            
            double rho = meas_package.raw_measurements_(0);
            double phi = meas_package.raw_measurements_(1);
            double rho_dot = meas_package.raw_measurements_(2);
            double x= rho * cos(phi);
            double y= rho * sin(phi);
            double vx = rho_dot * cos(phi);
            double vy = rho_dot * sin(phi);
            
            x_ << x, y, vx, vy;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            /**
             Initialize state.
             */
            float px = meas_package.raw_measurements_(0);
            float py = meas_package.raw_measurements_(1);
            
            x_(0) = px;
            x_(1) = py;
            
        }
        
        time_us_ = meas_package.timestamp_;
        // done initializing, no need to predict or update
        
        is_initialized_ = true;
        
        
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/
    
    
    //compute the time elapsed between the current and previous measurements
    float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;	//delta_t - expressed in seconds
    time_us_ = meas_package.timestamp_;
    
    if ( delta_t > 0.001 )
    {
        Prediction(delta_t);
    }
    

    /*****************************************************************************
     *  Update
     ****************************************************************************/
    
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        
        UpdateRadar(meas_package);
        
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
        
        UpdateLidar(meas_package);
    }
    
    }

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    
    /*****************************************************************************
     *  Create Augmented Sigma Points
     ****************************************************************************/
    
    //create augmented mean state
    
    x_aug_.head(5) = x_;
    x_aug_(5) = 0;
    x_aug_(6) = 0;
    
    //create augmented covariance matrix
    P_aug_.fill(0.0);
    P_aug_.topLeftCorner(5,5) = P_;
    P_aug_(5,5) = std_a_*std_a_;
    P_aug_(6,6) = std_yawdd_*std_yawdd_;
    
    //create square root matrix
    MatrixXd L = P_aug_.llt().matrixL();
    
    //create augmented sigma points
    Xsig_aug_.col(0)  = x_aug_;
    for (int i = 0; i< n_aug_; i++)
    {
        Xsig_aug_.col(i+1)       = x_aug_ + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * L.col(i);
    }
    
    
    /*****************************************************************************
     *  Predict Sigma Points
     ****************************************************************************/

    
    //predict sigma points
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
        //extract values for better readability
        double p_x = Xsig_aug_(0,i);
        double p_y = Xsig_aug_(1,i);
        double v = Xsig_aug_(2,i);
        double yaw = Xsig_aug_(3,i);
        double yawd = Xsig_aug_(4,i);
        double nu_a = Xsig_aug_(5,i);
        double nu_yawdd = Xsig_aug_(6,i);
        
        //predicted state values
        double px_p, py_p;
        
        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }
        
        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;
        
        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;
        
        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;
        
        //write predicted sigma point into right column
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
    }
    
        /*****************************************************************************
         *  Predict mean and covariance
         ****************************************************************************/

        
        //predicted state mean
        x_.fill(0.0);
        for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
            x_ = x_+ weights_(i) * Xsig_pred_.col(i);
        }
        
        //predicted state covariance matrix
        P_.fill(0.0);
        for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
            
            // state difference
            VectorXd x_diff = Xsig_pred_.col(i) - x_;
            //angle normalization
            x_diff(3) = NormalizeAngle(x_diff(3));
            
            P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
        }
    
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
    
    VectorXd z_pred = H_laser_ * x_;
    VectorXd y = meas_package.raw_measurements_ - z_pred;
    MatrixXd Ht = H_laser_.transpose();
    MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;
    
    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_laser_) * P_;
    
    lidar_NIS_ = y.transpose() * S.inverse()*y;
    
    cout<<"lidar_NIS: "<<lidar_NIS_<< endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
    
    
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        // extract values for better readibility
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v  = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;
        
        // measurement model
        if (p_x!=0 || p_y!=0){
            Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
            Zsig_(1,i) = atan2(p_y,p_x);                                 //phi
            Zsig_(1,i) = NormalizeAngle(Zsig_(1,i));
            Zsig_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);
        }//r_dot
    }
    
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z_);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig_.col(i);
    }
    
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z_,n_z_);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig_.col(i) - z_pred;
        
        //angle normalization
        z_diff(1) = NormalizeAngle(z_diff(1));
        
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    
    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z_,n_z_);
    R <<    std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0,std_radrd_*std_radrd_;
    S = S + R;

    
    /*****************************************************************************
     *  Update RADAR
     ****************************************************************************/

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z_);

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        //residual
        VectorXd z_diff = Zsig_.col(i) - z_pred;
        //angle normalization
        z_diff(1) = NormalizeAngle(z_diff(1));
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        x_diff(3) = NormalizeAngle(x_diff(3));
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    //residual
    VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
    
    //angle normalization
    z_diff(1) = NormalizeAngle(z_diff(1));
    
    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    
    radar_NIS_ = z_diff.transpose() * S.inverse()*z_diff;
    
    cout<<"radar_NIS : "<<radar_NIS_<< endl;

}

double UKF::NormalizeAngle(double angle){
    while (angle > M_PI) angle-=2.*M_PI;
    while (angle <-M_PI) angle+=2.*M_PI;
    return angle;
}


