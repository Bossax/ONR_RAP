DROP DATABASE IF EXISTS RAP;
CREATE DATABASE RAP;
USE RAP;
CREATE TABLE `HEM_original_position` (
  `tx_time` DATETIME NOT NULL,
  `longitude` decimal(11,8) NOT NULL,
  `latitude` decimal(11,8) NOT NULL,
  `altitude` decimal(5,3) NOT NULL,
  `heading` decimal(5,2) NOT NULL,
  `x_velocity` decimal(5,3) NOT NULL,
  `surface_range` decimal(6,3),
	`x_err` decimal(6,3) ,
	`y_err` decimal(6,3) ,
	`z_err` decimal(6,3) ,
	`traveltime_perturbation` decimal(5,3) NOT NULL,
	`SNR` decimal(4,2) NOT NULL,
    PRIMARY KEY (tx_time)
    );


