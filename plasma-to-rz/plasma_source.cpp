#include <iostream>
#include <vector>
#include <cmath>
#include "plasma_source.hpp"
#include <stdlib.h>     

std::vector<double> source_profile;
std::vector<double> ion_kt;


void setup_(const double *ion_density_ped, const double *ion_density_sep,
	    const double *ion_density_origin, const double *ion_temp_ped,
	    const double *ion_temp_sep, const double *ion_temp_origin, 
	    const double *pedistal_rad, const double *ion_density_peak,
	    const double *ion_temp_peak, const double *minor_radius, 
	    const double *major_radius, const double *elongation, 
	    const double *triangularity, const double *shafranov, 
	    const char* plasma_type, const int *plasma_id,
	    const int* number_of_bins ) {

  plasma plasma_parameters;

  //double random_store[4];

  // set params
  plasma_parameters.IonDensityPedistal = *ion_density_ped;
  plasma_parameters.IonDensitySeperatrix = *ion_density_sep;
  plasma_parameters.IonDensityOrigin = *ion_density_origin;
  plasma_parameters.IonTemperaturePedistal = *ion_temp_ped;
  plasma_parameters.IonTemperatureSeperatrix = *ion_temp_sep;
  plasma_parameters.IonTemperatureOrigin = *ion_temp_origin;
  plasma_parameters.PedistalRadius = *pedistal_rad;
  plasma_parameters.IonDensityPeaking = *ion_density_peak;
  plasma_parameters.IonTemperaturePeaking = *ion_temp_peak;
  plasma_parameters.MinorRadius = *minor_radius;
  plasma_parameters.MajorRadius = *major_radius;
  plasma_parameters.PlasmaType = plasma_type;
  plasma_parameters.PlasmaId = *plasma_id;

  setup_plasma_source(*number_of_bins,plasma_parameters);

  return;
}

/*
 * sample the pdf src_profile, to generate the sampled minor radius
 */
void sample_source_(const int *num_of_bins, const double *bin_width, double &sampled_radius, 
		    int &sampled_bin, double *rn_store1, double *rn_store2) {
  
  for ( int i = 0 ; i < *num_of_bins ; i++ ) {
    if ( *rn_store1 <= source_profile[i] ) {
      if ( i > 0 ) {
	sampled_radius = (float(i-1)*(*bin_width)) + (*bin_width*(*rn_store2));
	sampled_bin = i;
	return;
      } else {
	sampled_radius = *bin_width*(*rn_store2);
	sampled_bin = i;
	return;
      }
    }
  }

  std::cerr << "error" << std::endl;
  std::cerr << "Sample position greater than plasma radius" << std::endl;
  exit(1);
  return;
}

/*
 * sample the energy of the neutrons, updates energy neutron in mev
 */
void sample_energy_(const int *bin_number, double *random_number1, double *random_number2,
		    double &energy_neutron) {
  // generate the normally distributed number
  const double twopi = 6.28318530718;
  double sample1 = std::sqrt(-2.0*std::log(*random_number1));
  double sample2 = cos(twopi*(*random_number2));
  energy_neutron = (5.59/2.35)*(ion_kt[*bin_number])*sample1*sample2;
  energy_neutron += 14.08; 
  return;
}

/*
 * convert the sampled radius to an rz coordinate by using plasma parameters
 */
void convert_rad_to_rz_(const double *major_radius, const double *minor_radius,
			const double *elongation, const double *delta, 
			const double *shafranov, const double *minor_sampled,
			const double *rn_store, double &radius, double &height)
{
  const double twopi = 6.28318530718;
  
  double alpha = twopi*(*rn_store);
  //  std::cout << alpha << std::endl;
  double shift = *shafranov*(1.0-std::pow(*minor_sampled/(*minor_radius),2));
  //  std::cerr << std::pow(*minor_sampled/(*minor_radius),2) << std::endl;

  radius = (*major_radius) + (*minor_sampled)*cos(alpha+(*delta*sin(alpha))) + shift;
  height = (*elongation)*(*minor_sampled)*sin(alpha);
  
 
  return;
}
  

/*
 * convert rz_to_xyz
 */
void convert_rz_to_xyz_(const double *r, const double *rn_store, 
                        double &x, double &y, const double *min_toroidal_angle,
                        const double* max_toroidal_angle, double &toroidal_angle)
{
  double toroidal_extent = (*max_toroidal_angle)-(*min_toroidal_angle);
  toroidal_angle = toroidal_extent*(*rn_store)+(*min_toroidal_angle);
  x = *r*sin(toroidal_angle);
  y = *r*cos(toroidal_angle);
  return;
}

/*
 * sets up the cumulatitive probability profile
 * on the basis of the ion temp and ion density
 * this portion is deterministic
 */
void setup_plasma_source(int number_of_bins,
			 const plasma plasma_parameters )
{
  double ion_d; // ion density
  double ion_t; // ion temp
  double sig_dt; // dt xs

  std::vector<double> src_strength; // the source strength, n/m3
  double r;

  double bin_width = plasma_parameters.MinorRadius/float(number_of_bins);
  double total = 0.0; // total source strength

  for (int i = 0 ; i < number_of_bins ; i++) {
    r = bin_width * float(i);
    ion_d = ion_density(plasma_parameters,r);
    ion_t = ion_temperature(plasma_parameters,r);
    src_strength.push_back(std::pow(ion_d,2)*dt_xs(ion_t));
    ion_kt.push_back(sqrt(ion_t/1000.0)); // convert to sqrt(MeV)
    total += src_strength[i];
  }

  // normalise the source profile
  double sum =0 ;
  for ( int i = 0 ; i < number_of_bins ; i++) {
	sum += src_strength[i];
        source_profile.push_back(sum/total);
  }
  return;
}

/*
 * function that returns the ion density given the 
 * given the critical plasma parameters
 */
double ion_density(const plasma plasma_parameters,
		   double SampleRadius)
{
  double ion_dens = 0.0;

  if( plasma_parameters.PlasmaId == 0 ) {
    ion_dens = plasma_parameters.IonDensityOrigin*
      (1.0-std::pow(SampleRadius/plasma_parameters.MinorRadius,2));
  } else {
    if(SampleRadius <= plasma_parameters.PedistalRadius) {
      ion_dens += plasma_parameters.IonDensityPedistal;
      double product;
      product = 1.0-std::pow(SampleRadius/plasma_parameters.PedistalRadius,2);
      product = std::pow(product,plasma_parameters.IonDensityPeaking);
      ion_dens += (plasma_parameters.IonDensityOrigin-plasma_parameters.IonDensityPedistal)*
            	  (product);
    } else {
      ion_dens += plasma_parameters.IonDensitySeperatrix;
      double product;
      product = plasma_parameters.IonDensityPedistal - plasma_parameters.IonDensitySeperatrix;
      ion_dens += product*(plasma_parameters.MinorRadius-SampleRadius)/
	(plasma_parameters.MinorRadius-plasma_parameters.PedistalRadius);
    }
  }

  return ion_dens;
}
		   
/*
 * function that returns the ion density given the 
 * given the critical plasma parameters
 */
double ion_temperature(const plasma plasma_parameters,
		   double SampleRadius)
{
  double ion_temp = 0.0;

  if( plasma_parameters.PlasmaId == 0 ) {
    ion_temp = plasma_parameters.IonTemperatureOrigin*
      (1.0-std::pow(SampleRadius/plasma_parameters.MinorRadius,
		    plasma_parameters.IonTemperaturePeaking));
  } else {
    if(SampleRadius <= plasma_parameters.PedistalRadius) {
      ion_temp += plasma_parameters.IonTemperaturePedistal;
      double product;
      product = 1.0-std::pow(SampleRadius/plasma_parameters.PedistalRadius,2);
      product = std::pow(product,plasma_parameters.IonTemperaturePeaking);
      ion_temp += (plasma_parameters.IonTemperatureOrigin-
		   plasma_parameters.IonTemperaturePedistal)*(product);
    } else {
      ion_temp += plasma_parameters.IonTemperatureSeperatrix;
      double product;
      product = plasma_parameters.IonTemperaturePedistal - 
	plasma_parameters.IonTemperatureSeperatrix;
      ion_temp += product*(plasma_parameters.MinorRadius-SampleRadius)/
	(plasma_parameters.MinorRadius-plasma_parameters.PedistalRadius);
    }
  }

  return ion_temp;
}


/*
 * returns the dt cross section for a given ion temp
 */
double dt_xs(double ion_temp)
{
  double dt;
  double c[7]={2.5663271e-18,19.983026,2.5077133e-2,
	       2.5773408e-3,6.1880463e-5,6.6024089e-2,
	       8.1215505e-3};
  
  double u = 1.0-ion_temp*(c[2]+ion_temp*(c[3]-c[4]*ion_temp))
    /(1.0+ion_temp*(c[5]+c[6]*ion_temp));

  dt = c[0]/(std::pow(u,5./6.)*std::pow(ion_temp,2./3.0));
  dt *= exp(-1.*c[1]*std::pow(u/ion_temp,1./3.));

  return dt;
}
