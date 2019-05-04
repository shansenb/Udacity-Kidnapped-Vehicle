/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <map>

#include "helper_functions.h"

using std::string;
using std::vector;


void ParticleFilter::init(double x, double y, double theta, double std[]) {

  if (is_initialized) {
    return;
  }
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  std::default_random_engine gen;
  num_particles = 100;  // TODO: Set the number of particles
  std::normal_distribution<double> dist_x(x,std[0]);
  std::normal_distribution<double> dist_y(y,std[1]);
  std::normal_distribution<double> dist_theta(theta,std[2]);

  //initialize particles with normal distribution
  for (int i=0; i<num_particles; i++){
    Particle particle;
    particle.id=i;
    particle.x=dist_x(gen);
    particle.y=dist_y(gen);
    particle.theta=dist_theta(gen);
    particle.weight=1.0;

    particles.push_back(particle);
  }
  is_initialized=true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(0, std_pos[0]);
  std::normal_distribution<double> dist_y(0, std_pos[1]);
  std::normal_distribution<double> dist_theta(0, std_pos[2]);
  //initialize particles with normal distribution
  for (int i=0; i<num_particles; i++){
    double theta0=particles[i].theta;

    if (fabs(yaw_rate)>.000001){
      particles.at(i).x+=(velocity/yaw_rate)*(sin(theta0+yaw_rate*delta_t)-sin(theta0));
      particles.at(i).y+=(velocity/yaw_rate)*(-cos(theta0+yaw_rate*delta_t)+cos(theta0));
      particles.at(i).theta+=yaw_rate*delta_t;
    }
    else{
      particles.at(i).x+=velocity*cos(theta0)*delta_t;
      particles.at(i).y+=velocity*sin(theta0)*delta_t;
    }

    //add Gauss noise
    particles.at(i).x+=dist_x(gen);
    particles.at(i).y+=dist_y(gen);
    particles.at(i).theta+=dist_theta(gen);

  }
  
  
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  //nearest neighbor calc, O(mn)
  
  for (auto &obs: observations){
    double minObs=std::numeric_limits<double>::max();

    for (auto const &predict : predicted){
      double obsDist =dist(obs.x,obs.y,predict.x,predict.y);
      if (obsDist<minObs){
        minObs=obsDist;
        obs.id=predict.id;

      }
      
    }

  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  //based on particle location, predict location of all map landmarks within sensor range


  double sig_x = std_landmark[0];
  double sig_y = std_landmark[1];
  double gauss_norm= 1 / (2 * M_PI * sig_x * sig_y);
  
  for (auto &particle: particles){


    vector<LandmarkObs> predicted;

    //for each landmark 
    for(auto const &landmark: map_landmarks.landmark_list){
      if (dist(landmark.x_f,landmark.y_f,particle.x,particle.y)<sensor_range){
        //we have a landmark that may be measured
        predicted.push_back(LandmarkObs{landmark.id_i,landmark.x_f,landmark.y_f});
      }
    }

    vector<LandmarkObs> localObs;
    
    //transform observations to global coords
    for (const auto &obs: observations){
      double ox=obs.x*cos(particle.theta)-obs.y*sin(particle.theta)+particle.x;
      double oy=obs.x*sin(particle.theta)+obs.y*cos(particle.theta)+particle.y;
      localObs.push_back({obs.id,ox,oy});
    }

    dataAssociation(predicted,localObs);


    //update the weight 
    particle.weight=1.0;

    for (auto const &obs: localObs){
      //update weight for particle based on conditional probability
      for( auto const &predict: predicted){
        if (predict.id==obs.id){
          
            // calculate exponent
            double exponent = (pow(obs.x - predict.x, 2) / (2 * pow(sig_x, 2)))
                        + (pow(obs.y - predict.y, 2) / (2 * pow(sig_y, 2)));
            
            // calculate weight using normalization terms and exponent
            particle.weight *= gauss_norm * exp(-exponent);
        }
      }        
    }



  }
  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  std::default_random_engine gen;
  vector<Particle> resampled;
  vector<double> weights;
  for (auto particle: particles) {
    weights.push_back(particle.weight);
  }

  // generate random starting index for resampling wheel
  std::uniform_int_distribution<int> uniformIntDist(0, num_particles-1);
  auto idx = uniformIntDist(gen);
  
  double max_weight = *(max_element(weights.begin(), weights.end()));
  std::uniform_real_distribution<double> uniformRealDist(0.0, max_weight);

  double beta = 0.0;
  for (int i = 0; i < num_particles; i++) {
    beta += 2*uniformRealDist(gen);
    while (beta>weights[idx]) {
      beta -= weights[idx];
      idx = (idx + 1) % num_particles;
    }
    resampled.push_back(particles[idx]);
  }

  particles = resampled;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}