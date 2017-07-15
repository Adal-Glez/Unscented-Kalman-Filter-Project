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
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;
    
  //define spreading parameter
  lambda_ = 3 - n_aug_;
  
  n_sigma_points =2*n_aug_+1;
    
  ///* augmented sigma points matrix
  Xsig_aug_ = MatrixXd(n_aug_, n_sigma_points);
   
  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_points);
    
  /// set vector for weights
    
  weights_ = VectorXd(n_sigma_points);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;

  for (int i=1; i<n_sigma_points; i++) {
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
   
    float delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
    
    previous_timestamp_ = meas_package.timestamp_;
    
    Prediction(delta_t);
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        UpdateRadar(meas_package);
    } else {
        // Laser
        UpdateLidar(meas_package);
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
    AugmentedSigmaPoints();
    SigmaPointPrediction(delta_t);
    PredictMeanAndCovariance();
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
    int n_z = meas_package.raw_measurements_.rows();
    
    VectorXd z_pred(n_z);
    MatrixXd S(n_z, n_z);
    
    //MatrixXd Xsig_aug(n_aug_, n_aug_sigma_);
    PredictLidarMeasurement();
    UpdateState(meas_package);
    
    // calculate NIS
    VectorXd z_diff = meas_package.raw_measurements_-z_pred;
    NIS_lidar_ = z_diff.transpose() * S.inverse() * z_diff;

}

void UKF::PredictLidarMeasurement() {
    
    //transform sigma points into measurement space
    for(int i=0; i < (2*n_aug_ + 1); i++){
        double px = Xsig_pred_.col(i)(0);
        double py = Xsig_pred_.col(i)(1);
        
        Zsig_lidar_.col(i) << px, py;
        
    }
    
    //calculate mean predicted measurement
    z_pred_lidar_.fill(0.);
    for(int i=0; i < (2*n_aug_ + 1); i++){
        float weight = weights_(i);
        z_pred_lidar_ += weight * Zsig_lidar_.col(i);
    }
    
    S.fill(0);
    //calculate measurement covariance matrix S
    for(int i=0; i < (n_sigma_points); i++){
        float weight = weights_(i);
        VectorXd z_diff = Zsig_lidar_.col(i) - z_pred_lidar_;
        S_lidar_ += weight * z_diff * z_diff.transpose();
    }
    
    MatrixXd R = MatrixXd(n_z_lidar_, n_z_lidar_);
    R.fill(0);
    R(0,0) = std_laspx_ * std_laspx_;
    R(1,1) = std_laspy_ * std_laspy_;
    
    S = S + R;

}

void UKF::UpdateState(MeasurementPackage meas_package) {
    
    //create matrix for cross correlation Tc
    int n_z = z_pred.rows();
    int n_x = x.rows();
    
    MatrixXd Tc = MatrixXd(n_x, n_z);
    
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < n_aug_sigma_; i++) {  //2n+1 simga points
        
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        // state difference
        VectorXd x_diff = Xsig_pred.col(i) - x;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    //residual
    VectorXd z_diff = z - z_pred;
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    //update state mean and covariance matrix
    x = x + K * z_diff;
    P = P - K*S*K.transpose();
    
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
    int n_z = meas_package.raw_measurements_.rows();
    
    VectorXd z_pred(n_z);
    MatrixXd S(n_z, n_z);
    
    PredictRadarMeasurement();
    UpdateState(meas_package);
    
    VectorXd z_diff = meas_package.raw_measurements_-z_pred;
    
    while (z_diff(2)> M_PI) z_diff(2)-=2.*M_PI;
    while (z_diff(2)<-M_PI) z_diff(2)+=2.*M_PI;
    
    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}

void UKF::PredictRadarMeasurement() {
    //  //mean predicted measurement
    //  VectorXd z_pred = VectorXd(n_z_radar_);
    //
    //  //measurement covariance matrix S
    //  MatrixXd S = MatrixXd(n_z_radar_,n_z_radar_);
    
    //transform sigma points into measurement space
    for(int i=0; i < (2*n_aug_ + 1); i++){
        double px = Xsig_pred_.col(i)(0);
        double py = Xsig_pred_.col(i)(1);
        double v  = Xsig_pred_.col(i)(2);
        double psi = Xsig_pred_.col(i)(3);
        
        double rou = sqrt(px*px + py*py);
        double phi = atan2(py, px);
        double rou_dot;
        
        if(rou < 0.0000001){
            rou_dot = 0.0;
        }else{
            rou_dot = 1/rou * (px*cos(psi)*v + py*sin(psi)*v);
        }
        
        Zsig_radar_.col(i) << rou , phi , rou_dot;
        
    }
    
    //calculate mean predicted measurement
    z_pred_radar_.fill(0.);
    for(int i=0; i < (n_sigma_points); i++){
        float weight = weights_(i);
        z_pred_radar_ += weight * Zsig_radar_.col(i);
    }
    
    S.fill(0);
    //calculate measurement covariance matrix S
    for(int i=0; i < (2*n_aug_ + 1); i++){
        float weight = weights_(i);
        VectorXd z_diff = Zsig_radar_.col(i) - z_pred_radar_;
        S += weight * z_diff * z_diff.transpose();
    }
    
    MatrixXd R = MatrixXd(n_z_radar_, n_z_radar_);
    R.fill(0);
    R(0,0) = std_radr_ * std_radr_;
    R(1,1) = std_radphi_ * std_radphi_;
    R(2,2) = std_radrd_ * std_radrd_;
    
    S = S+R;
    

}

void UKF::AugmentedSigmaPoints() {
    
    //create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);
    
    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(7, 7);
    
    
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
    for (int i = 0; i < n_aug_; i++)
    {
        Xsig_aug.col(i+1)     = x_aug + sqrt(lambda_ +n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda_ +n_aug_) * L.col(i);
    }
   
}

void UKF::SigmaPointPrediction(float delta_t) {
    
    //predict sigma points
    for (int i = 0; i< n_sigma_points; i++)
    {
        //extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);
        
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
        Xsig_pred(0,i) = px_p;
        Xsig_pred(1,i) = py_p;
        Xsig_pred(2,i) = v_p;
        Xsig_pred(3,i) = yaw_p;
        Xsig_pred(4,i) = yawd_p;
    }
    
    
}

void UKF::PredictMeanAndCovariance() {
    
    //predicted state mean
    x.fill(0.0);
    for (int i = 0; i < n_sigma_points; i++) {  //iterate over sigma points
        x = x+ weights(i) * Xsig_pred.col(i);
    }
    
    //predicted state covariance matrix
    P.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        
        // state difference
        VectorXd x_diff = Xsig_pred.col(i) - x;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        P = P + weights(i) * x_diff * x_diff.transpose() ;
    }
    
}



