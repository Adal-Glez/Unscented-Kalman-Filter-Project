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

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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
  //set state dimension
  n_x = 5;

  //set augmented dimension
  n_aug = 7;
    
  //define spreading parameter
  lambda = 3 - n_aug;
  
  n_sigma_points =2*n_aug+1;
    
    ///* augmented sigma points matrix
    Xsig_aug_ = MatrixXd(n_aug_, n_sigma_points);
   
    ///* predicted sigma points matrix
    Xsig_pred_ = MatrixXd(n_x_, n_sigma_points);
    
    /// set vector for weights
    
    weights_ = VectorXd(n_sigma_points_);
    double weight_0 = lambda_ / (lambda_ + n_aug);
    weights_(0) = weight_0;

    for (int i=1; i<n_sigma_points; i++) {
        double weight = 0.5/(n_aug+lambda);
        weights_(i) = weight;
    }
    
    

    
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
    if (!is_initialized_) {
        
        //x_ << 1, 1, 1, 1;
        P_ << 1., 0., 0., 0.,0.,
            0., 1., 0., 0., 0.,
            0., 0., 1., 0., 0.,
            0., 0., 0., 1., 0.,
            0., 0., 0., 0., 1.;
        
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            
            // convert from polar to Cartesian
            float rho = meas_package.raw_measurements_[0];
            float phi = meas_package.raw_measurements_[1];
            
            x_ << rho * cos(phi), rho * sin(phi), 0, 0, 0;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            
            x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
        }
        // initialize the timestamp
        previous_timestamp_ = meas_package.timestamp_;
        
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;

    }
   
    float dt = (meas_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
    
    previous_timestamp_ = meas_pack.timestamp_;
    
    Prediction(dt);
    
    if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        UpdateRadar(meas_package);
    } else {
        // Laser
        ekf_.Update(meas_package);
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
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
    MatrixXd Xsig_aug(n_aug_, n_aug_sigma_);
    AugmentedSigmaPoints(x_, P_, Xsig_aug);
    sigmapoints prediction
    predict mean and cov
}


void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
    
    //set state dimension
    //int n_x = 7;
    
    //set augmented dimension
    int n_aug = 7;
    
    //Process noise standard deviation longitudinal acceleration in m/s^2
    double std_a = 0.2;
    
    //Process noise standard deviation yaw acceleration in rad/s^2
    double std_yawdd = 0.2;
    
    //define spreading parameter
    double lambda = 3 - n_aug;
    
    //create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug);
    
    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(7, 7);
    
    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
    
    
    //create augmented mean state
    x_aug.head(5)=x;
    x_aug(5)=0;
    x_aug(6)=0;
  
    //create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5,5) = P;
    P_aug(5,5)= std_a * std_a;
    P_aug(6,6)=std_yawdd * std_yawdd;
    
    
    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();
    
    //create augmented sigma points
    Xsig_aug.col(0) = x_aug;
    for (int i = 0; i < n_aug; i++)
    {
        Xsig_aug.col(i+1)     = x_aug + sqrt(lambda+n_aug) * L.col(i);
        Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda+n_aug) * L.col(i);
    }
   


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
}
