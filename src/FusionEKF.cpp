#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  
  //Lesson 24:11
  H_laser_ << 1, 0, 0, 0,
  			 0, 1, 0, 0;

  /**
   * From Lesson 24:8
   * DONE: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  
  //object state
  VectorXd x = VectorXd(4);
  x << 0, 0, 0, 0;
  
  //object covariance matrix
  MatrixXd P = MatrixXd(4,4);
  P << 0, 0, 0, 0,
       0, 0, 0, 0,
       0, 0, 0, 0,
       0, 0, 0, 0;
  
  //external motion - I'm not sure if we need this. 
  VectorXd u = VectorXd(4);
  u << 0, 0, 0, 0;
  
  
  //state transition matrix Lesson 24:13
  MatrixXd F = MatrixXd(4,4);
  F << 1, 0, 0.5, 0,
       0, 1, 0, 0.5,
       0, 0, 1, 0,
       0, 0, 0, 1;
  
  //process covariance matrix
  MatrixXd Q = MatrixXd(4,4);
  Q << 0, 0, 0, 0,
       0, 0, 0, 0,
       0, 0, 0, 0,
       0, 0, 0, 0;
  
  /**
   * Init Initializes Kalman filter
   * @param x_in Initial state
   * @param P_in Initial state covariance
   * @param F_in Transition matrix
   * @param H_in Measurement matrix
   * @param R_in Measurement covariance matrix
   * @param Q_in Process covariance matrix
   */
  
  ekf_.Init( x, P, F, H_laser_, R_laser_, Q);
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * DONE: Initialize the state ekf_.x_ with the first measurement.
     * DONE: Create the covariance matrix.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      
      //Convert radar from polar to cartesian coordinates 
      // 0: r, 1: theta, 2: r dot
      float px = cos(measurement_pack.raw_measurements_[1]) * measurement_pack.raw_measurements_[0];
      float py = sin(measurement_pack.raw_measurements_[1]) * measurement_pack.raw_measurements_[0];
      float vx = cos(measurement_pack.raw_measurements_[1]) * measurement_pack.raw_measurements_[2]; 
      float vy = sin(measurement_pack.raw_measurements_[1]) * measurement_pack.raw_measurements_[2];
      
      //and initialize state covarience
      ekf_.P_ << 1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 1, 0,
                 0, 0, 0, 1;
      
      //and initialize state
      ekf_.x_<< px, py, vx, vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state.
    
      float px = measurement_pack.raw_measurements_[0];
      float py = measurement_pack.raw_measurements_[1];
      
      //Lesson 24:14 tracking.cpp
      ekf_.P_ << 1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 1000, 0,
                 0, 0, 0, 1000;
      
      ekf_.x_ << px, py, 0,0;
    }
	
    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }

  /**
   * Prediction
   */
	
  /**
   * TODO: 
   * 
   */
 
  //Lesson 24:13
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //seconds
  previous_timestamp_ = measurement_pack.timestamp_;
 
  //Update the state transition matrix F according to the new elapsed time.
  //update the dt values in F
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;
  
  //Update the process noise covariance matrix. Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
  
  float noise_ax = 9;
  float noise_ay = 9;
  
  float dtTo2 = dt * dt;
  float dtTo3 = dt * dt * dt;
  float dtTo4 = dt * dt * dt * dt;
  
  //From Lesson 24:10
  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << dtTo4/4 * noise_ax, 0, dtTo3/2 * noise_ax, 0,
  			0, dtTo4/4 * noise_ay, 0, dtTo3/2 * noise_ay,
  			dtTo3/2 * noise_ax, 0, dtTo2 * noise_ax, 0,
  			0, dtTo3/2 * noise_ay, 0, dtTo2 * noise_ay;
  
  
  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    
  } else {
    // Laser updates
	ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

/*
You will need to:
	1. initialize variables and matrices (x, F, H_laser, H_jacobian, P, etc.)
	2. initialize the Kalman filter position vector with the first sensor measurements
	3. modify the F and Q matrices prior to the prediction step based on the elapsed time between measurements
	4. call the update step for either the lidar or radar sensor measurement. Because the update step for 
            lidar and radar are slightly different, there are different functions for updating lidar and radar.
*/
