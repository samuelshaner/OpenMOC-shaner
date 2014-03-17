/**
 * @file FunctionalMaterial.h
 * @brief
 * @date September 14, 2013
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef FUNCTIONALMATERIAL_H_
#define FUNCTIONALMATERIAL_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <math.h>
#include "Material.h"
#include <vector>
#include "TimeStepper.h"
#endif


/**
 * @class FunctionalMaterial FunctionalMaterial.h "openmoc/src/host/FunctionalMaterial.h"
 * @brief The material class represents a unique material and its relevant
 *        nuclear data (ie, multigroup cross-sections) for neutron transport.
 *        The cross sections are allowed to have time or temperature dependence.
 */
class FunctionalMaterial : public Material {

protected:

    /** An array of the precursor concentrations */
    double* _prec_conc;

    /** An array of the precursor frequencies */
    double* _prec_freq;

    /** The material temperature */
    double* _temperature;

    /** The time grid used for interpolating cross sections */
    double* _time;

    /** The number of delayed neutron energy groups */
    int _num_delay_groups;

    double* _gamma;

    bool _sigma_a_func_temp;
    bool _sigma_a_func_time;
    bool _sigma_s_func_time;

    /* flag to adjust inscatter to conserve total xs */
    bool _conserve_sigma_t;

    int _num_time_steps;

    TimeStepper* _ts;

public:
    FunctionalMaterial(short int id);
    virtual ~FunctionalMaterial();

    /* set number of energy groups */
    virtual void setNumEnergyGroups(const int num_groups, const int num_time_steps);

    /* set sigma a */
    void setSigmaATime(int num_time_steps, int num_groups, double* xs);

    /* set sigma s */
    void setSigmaSTime(int num_time_steps, int num_groups, double* xs);
    
    /* copy and clone material */
    virtual FunctionalMaterial* clone();

    /* get and set time */
    void setTime(double* time, int num_time_steps);
    double* getTime();

    /* set bool for temp functional properties */
    void sigmaAFuncTemp(bool func_temp);

    /* set bool for time functional properties */
    virtual void setSigmaA(double* xs, int num_groups);
    virtual void setSigmaS(double* xs, int num_groups);
    void sigmaAFuncTime(bool func_time);
    void sigmaSFuncTime(bool func_time);

    /* set and get gamma, the doppler feedback coefficient */
    void setGamma(double* gamma, int num_groups);
    double* getGamma();

    /* sync material properties to a particular state */
    void sync(materialState state);


    void copySigmaSRef(Material* material);

    /* interpolate a particular cross section */
    double interpolateXS(double* xs_vec, materialState state, int group);
    double interpolateScatterXS(double* xs_vec, materialState state, int group_from, int group_to);

    /* get and set transient material properties */
    void setPrecConc(materialState state, double conc, int group);
    void setPrecFreq(materialState state, double conc, int group);
    double getPrecConc(materialState state, int group);
    double getPrecFreq(materialState state, int group);
    void copyPrecConc(materialState state_from, materialState state_to);
    void copyPrecFreq(materialState state_from, materialState state_to);
    void initializeTransientProps(double num_delay_groups, bool cmfd_mesh);

    /* set time stepper */
    void setTimeStepper(TimeStepper* ts);

    void setConserveSigmaT(bool conserve_sigma_t);
};


#endif /* FUNCTIONALMATERIAL_H_ */
