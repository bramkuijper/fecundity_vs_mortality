// Gillespie simulation 
// for the evolution of parental effects
// given mortality vs fecundity selection

#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "auxiliary.h"

//#define NDEBUG

using namespace std;


// random number generator 
// see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

const size_t Npatch = 200; 
const size_t Npp = 100;

unsigned long int NumGen = 5e08; // maximum number of generations
int sample = 20;
double d0 = 0.01; // starting value dispersal
double h0 = 0.01; // starting value phenotype switching

// initial values for the parental effect loci
double p1_init = 0.0;
double p2_init = 0.0;

double error_prob = 0.5; // probability that phenotype is something else...

int seed = -1; // the random seed
time_t total_time;  // track the total time the simulation is running


unsigned long int generation = 0; // track the current generation 

// output once every x generations
unsigned long int output_interval = 2e04;

double envt_switch[2] = {0,0}; // patch state switch rate
double k = 0.001; // cost of dispersal

// mortality of maladapted breeders in both envts
double c[2] = {0,0}; 

// mortality of adapted breeders in both envts
double C[2] = {0,0}; 

// fecundity of maladapted breeders in both envts
double f[2] = {0,0}; 

// fecundity of adapted breeders in both envts
double F[2] = {0,0}; 

// mutation rates
double mu_h = 0.01;
double sdmu_h = 0.01;

// diploid individuals
struct Individual
{
    double d[2]; // one locus coding for migration
    double p1[2]; // proportion z1 offspring by z1 mother 
    double p2[2]; // proportion z1 offspring by z2 mother
    bool z; // state of the individual wrt its environment, 0, 1
};

struct Patch
{
    Individual locals[Npp]; // the local senior
    bool state; //environmental state of the patch

};

Patch MetaPop[Npatch];


// statistics for sampling
// number of patches in environment one or two
// with 0, 1 to Npp maladapted individuals
// hence a total of Npp+1 combinations
size_t Npatches[2][Npp+1];

// next to counts store lists with ids
size_t patch_ids[2][Npp+1][Npatch];

// give the outputfile a unique name
string filename("sim_px");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

// initialize the command line arguments to vary 
// the parameters
void init_arguments(int argc, char *argv[])
{
    envt_switch[0] = atof(argv[1]);
    envt_switch[1] = atof(argv[2]);
    c[0] = atof(argv[3]);
    c[1] = atof(argv[4]);
    C[0] = atof(argv[5]);
    C[1] = atof(argv[6]);
    f[0] = atof(argv[7]);
    f[1] = atof(argv[8]);
    F[0] = atof(argv[9]);
    F[1] = atof(argv[10]);
    k = atof(argv[11]);
    mu_h = atof(argv[12]);
    sdmu_h = atof(argv[13]);
    d0 = atof(argv[14]);
    h0 = atof(argv[15]);
    error_prob = atof(argv[16]);
    p1_init = atof(argv[17]);
    p2_init = atof(argv[18]);
}

// initialization function, runs at start of sim
// gives all the individuals some starting
// value
void init_pop()
{ 
    // start the time
    total_time = time(NULL);

    // obtain a seed from current nanosecond count
	seed = get_nanoseconds();

    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);

    // initialize patch state counters and set them to 0
    for (int envt_i = 0; envt_i < 2; ++envt_i)
    {
        for (int n_adapted = 0; n_adapted <= Npp; ++n_adapted)
        {
            Npatches[envt_i][n_adapted] = 0;
        }
    }

    size_t n_adapted;

    // initialize all the patches
    for (size_t i = 0; i < Npatch; ++i)
    {
        //patches are randomly put in env + or -
        MetaPop[i].state = gsl_rng_uniform(r) < 0.5; 
    
        n_adapted = 0;

        //initialize breeders
        for (size_t j = 0; j < Npp; ++j)
        {
            // randomly assign initial phenotypes
            MetaPop[i].locals[j].z = gsl_rng_uniform(r) < 0.5;

            if (MetaPop[i].locals[j].z == MetaPop[i].state)
            {
                ++n_adapted;
            }
        
            // initialize allelic values
            for (int allele = 0; allele < 2; ++allele)
            {
                MetaPop[i].locals[j].d[allele] = d0;
                MetaPop[i].locals[j].p1[allele] = p1_init;
                MetaPop[i].locals[j].p2[allele] = p2_init;
            }
        }

        // increase the counter of patches of
        // (i) environmental state 'state'
        // (iii) containing n_adapted individuals
        //
        size_t Npatches_of_this_type = 
            Npatches[MetaPop[i].state][n_adapted];

        // append this patch to the corresponding stack of patch id's
        patch_ids[MetaPop[i].state][n_adapted][Npatches_of_this_type] = i;

        // update counter
        ++Npatches[MetaPop[i].state][n_adapted];
    }
}


// mutate allelic values
double mut(double p)
{
    assert(p >= 0 && p <= 1);
    p+= gsl_rng_uniform(r) < mu_h ? gsl_ran_gaussian(r,sdmu_h) : 0;
    p = p < 0 ? 0 : p > 1 ? 1 : p;

    return(p);
}

// make a kid out of parent asexually
// add it either to local offspring pool
// or to pool of dispersers
void create_kid(Individual &mother, Individual &Kid)
{
    // if it is a phenotype 1 mother, 
    // we'll need to express p1
    // otherwise p2
    double expr_h;

    bool mom_has_z1 = mother.z == 0;

    if (gsl_rng_uniform(r) < error_prob)
    {
        mom_has_z1 = !mom_has_z1;
    }

    // mother's phenotype is one
    expr_h = mom_has_z1 ? 0.5 * (mother.p1[0] + mother.p1[1]) 
                                :
                                0.5 * (mother.p2[0] + mother.p2[1]);

    // offspring gets z1 or z2 phenotype
    Kid.z = gsl_rng_uniform(r) < expr_h ? 0 : 1;

//  ok such an error does not matter too much
//    if (gsl_rng_uniform(r) < error_prob)
//    {
//        Kid.z = !Kid.z;
//    }

    // inherit alleles
    for (int allele_i = 0; allele_i < 2; ++allele_i)
    {
        // turned of mutation on d
        Kid.d[allele_i] = mother.d[allele_i];

        // inherit p1 and p2
        Kid.p1[allele_i] = mut(mother.p1[allele_i]);
        Kid.p2[allele_i] = mut(mother.p2[allele_i]);
    }
}

// do the stats on the data
void write_data()
{
    // get variance and means
    double p1,p2,d,vard,varp1,varp2;
    double meanp1 = 0;
    double meanp2 = 0;
    double ssp1 = 0;
    double ssp2 = 0;

    double meand = 0;
    double ssd = 0;

    // loop through patches and individuals
    // and get stats on genotypes and patch freqs.
    for (size_t i = 0; i < Npatch; ++i)
    {
        for (size_t j = 0; j < Npp; ++j)
        {
            d = 0.5 * (MetaPop[i].locals[j].d[0] + MetaPop[i].locals[j].d[1]);
            p1 = 0.5 * (MetaPop[i].locals[j].p1[0] + MetaPop[i].locals[j].p1[1]);
            p2 = .5 * (MetaPop[i].locals[j].p2[0] + MetaPop[i].locals[j].p2[1]);

            meanp1+=p1;
            meanp2+=p2;
            meand+=d;
            ssp1+=p1*p1;
            ssp2+=p2*p2;
            ssd+=d*d;
        }    
    }

    DataFile << generation << ";"; 

    // write down all the patch frequency combinations
    for (size_t envt_i = 0; envt_i < 2; ++envt_i)
    {
        for (size_t n_adapted = 0; n_adapted <= Npp; ++n_adapted)
        {
            DataFile << setprecision(5) << (double)Npatches[envt_i][n_adapted] / Npatch << ";";
        }
    }

    meanp1 /= Npatch * Npp;
    meanp2 /= Npatch * Npp;
    meand /= Npatch * Npp;
    vard = ssd / (Npatch * Npp) - meand*meand;
    varp1 = ssp1 / (Npatch * Npp) - meanp1*meanp1;
    varp2 = ssp2 / (Npatch * Npp) - meanp2*meanp2;

    DataFile 
        << setprecision(5) << meanp1 << ";"
        << setprecision(5) << meanp2 << ";"
        << setprecision(5) << varp1 << ";"
        << setprecision(5) << varp2 << ";"
        << setprecision(5) << meand << ";"
        << setprecision(5) << vard << endl;
}

// headers of the datafile
void write_data_headers()
{
    DataFile << "generation;";

    for (size_t envt_i = 0; envt_i < 2; ++envt_i)
    {
        for (size_t n_adapted = 0; n_adapted <= Npp; ++n_adapted)
        {
            DataFile << setprecision(5) << "f_" << (envt_i+1) << "_" << (n_adapted) << ";";
        }
    }
    
    DataFile << "meanp1;meanp2;varp1;varp2;meand;vard" << endl;
}

// parameters at the end of the sim
void write_parameters()
{
    total_time = time(NULL) - total_time;
    DataFile << endl << endl << "d0;" << d0 << endl
                << "type;" << "phenotype-dependent" << endl
                << "p0;" << h0 << endl
                << "mu_h;" << mu_h << endl
                << "sdmu_h;" << sdmu_h << endl
                << "error_prob;" << error_prob << endl
                << "numgen;" << NumGen << endl
                << "k;" << k << endl
                << "s1;" << envt_switch[0] << endl
                << "s2;" << envt_switch[1] << endl
                << "c1;" << c[0] << endl
                << "c2;" << c[1] << endl
                << "C1;" << C[0] << endl
                << "C2;" << C[1] << endl
                << "f1;" << f[0] << endl
                << "f2;" << f[1] << endl
                << "F1;" << F[0] << endl
                << "F2;" << F[1] << endl
                << "seed;" << seed << endl
                << "p1_init;" << p1_init << endl
                << "p2_init;" << p2_init << endl
                << "Npatch;" << Npatch << endl
                << "Npp;" << Npp << endl
                << "runtime;" << total_time << endl;
}

void switch_patch_state(size_t patch_envt, size_t n_adapted)
{
    // sample a random patch that fulfills the conditions
    size_t random_patch = 
        gsl_rng_uniform_int(r, 
            Npatches[patch_envt][n_adapted]
        );

    // get the patch_id 
    size_t random_patch_id = 
        patch_ids[patch_envt][n_adapted][random_patch];

    // error checking
    assert(MetaPop[random_patch_id].state == patch_envt);
   
#ifndef NDEBUG
    size_t n_adapted_check = 0;

    for (size_t breeder_i = 0; breeder_i < Npp; ++breeder_i)
    {
        n_adapted_check += MetaPop[random_patch_id].locals[breeder_i].z 
            == MetaPop[random_patch_id].state;
    }

    assert(n_adapted_check == n_adapted);
#endif

    // update patch characteristics
    bool patch_envt_new = !MetaPop[random_patch_id].state;
    MetaPop[random_patch_id].state = patch_envt_new;
        
    size_t n_adapted_new = Npp - n_adapted; 

    // update statistics
    //
    // delete patch id in the corresponding stack
    // by replacing it with the last patch id in the stack
    // and then deleting the last element (by reducing the counter by 1)
    patch_ids[patch_envt][n_adapted][random_patch] =
        patch_ids[patch_envt][n_adapted][
            --Npatches[patch_envt][n_adapted] 
        ];

    // add patch id to the correct stack
    patch_ids[patch_envt_new][n_adapted_new][
        Npatches[patch_envt_new][n_adapted_new]++
    ] = random_patch_id;
    
    // error checking
    assert(MetaPop[random_patch_id].state == !patch_envt);
   
#ifndef NDEBUG
    n_adapted_check = 0;

    for (size_t breeder_i = 0; breeder_i < Npp; ++breeder_i)
    {
        n_adapted_check += MetaPop[random_patch_id].locals[breeder_i].z 
            == MetaPop[random_patch_id].state;
    }

    assert(n_adapted_check == n_adapted_new);
#endif
}

// mortality of an individual in one of the patches in state 
// patch envt
//
// size_t patch_envt: the environment of the patch the mortality event took place in
// size_t n_adapted: a patch containing n_adapted adapted individuals
// bool mortality_maladapted: true, maladapted guy dies; false, adapted guy dies;

void mortality(size_t patch_envt, size_t n_adapted, bool mortality_maladapted)
{
    // sample a random patch that fulfills the conditions
    size_t random_patch = 
        gsl_rng_uniform_int(r, 
            Npatches[patch_envt][n_adapted]
        );

    // get the patch_id 
    size_t random_patch_id = 
        patch_ids[patch_envt][n_adapted][random_patch];

    // error checking
    assert(MetaPop[random_patch_id].state == patch_envt);
   
#ifndef NDEBUG
    size_t n_adapted_check = 0;

    for (size_t breeder_i = 0; breeder_i < Npp; ++breeder_i)
    {
        n_adapted_check += MetaPop[random_patch_id].locals[breeder_i].z 
            == MetaPop[random_patch_id].state;
    }

    assert(n_adapted_check == n_adapted);
#endif

    // fecundity selection absent
    // so whether we get immigrant 
    // only depends on dispersal probability
    
    size_t sampled = Npp + sample;

    // array to store in cumulative distribution over all sampled
    // individuals
    double choose_prob[sampled];

    // array to store patches from which the sampled individuals
    // come from
    size_t patch_origin[sampled];

    // array to store the ids of the individuals
    size_t individual_ids[sampled];

    // array with the individuals that may die
    size_t mortality_candidate_ids[Npp];
    size_t n_mortality_candidates = 0;


    double sumprob = 0;
    double d = 0;
    double fecundity = 0;

    // make cumulative distribution of local competitive
    // ability over the locals (as dispersal is the only
    // trait affecting local competition just make a cumul
    // distribution of dispersal probabilities
    //
    // at the same time calculate which local breeders
    // might die
    for (size_t local_i = 0; local_i < Npp; ++local_i)
    {
        // calculate a local offspring's dispersal
        // capability
        d = 0.5 * (
                MetaPop[random_patch_id].locals[local_i].d[0] + 
                MetaPop[random_patch_id].locals[local_i].d[1]
                );

        // store the data
        patch_origin[local_i] = random_patch_id;
        individual_ids[local_i] = local_i;

        // calculate fecundity payoff
        // default: individual maladapted
        fecundity = f[patch_envt];

        // update fecundity payoff if adapted
        if (MetaPop[random_patch_id].locals[local_i].z == patch_envt)
        {
            fecundity = f[patch_envt];
        }

        // store the cumulative distribution
        choose_prob[local_i] = sumprob + (1-d) * fecundity;
        sumprob = choose_prob[local_i];

        // store the mortality candidates
        if (
                (mortality_maladapted && 
                    MetaPop[random_patch_id].locals[local_i].z 
                        != patch_envt
                )
                ||
                (!mortality_maladapted && 
                    MetaPop[random_patch_id].locals[local_i].z 
                        == patch_envt
                )

            )
        {
            mortality_candidate_ids[n_mortality_candidates++] = local_i;
        }
    }

    // we need to make sure that a very large sample of dispersers
    // receives equal sampling weighting relative to the Npp local individuals
    double pi = (double) Npp / sample;

    size_t id = 0;

    // make cumulative distribution over all sampled remote offspring
    // mimicking competition for breeding spots
    for (size_t remote_i = Npp; remote_i < sampled; ++remote_i)
    {
        // randomly pick patch
        patch_origin[remote_i] = gsl_rng_uniform_int(r,Npatch);

        // randomly pick individual within patch
        id = gsl_rng_uniform_int(r,Npp);

        individual_ids[remote_i] = id;

        // store individual id
        individual_ids[remote_i] = id;

        fecundity = f[MetaPop[patch_origin[remote_i]].state];

        if (MetaPop[patch_origin[remote_i]].locals[id].z ==
             MetaPop[patch_origin[remote_i]].state)
        {
            fecundity = F[MetaPop[patch_origin[remote_i]].state];
        }

        d =  0.5 * (
                        MetaPop[patch_origin[remote_i]].locals[id].d[0] + 
                        MetaPop[patch_origin[remote_i]].locals[id].d[1]
                        );

        choose_prob[remote_i] = sumprob + (1-k) * pi * d * fecundity;
        sumprob = choose_prob[remote_i];
    }

    double deviate = gsl_rng_uniform(r) * sumprob;

    Individual Kid;

    size_t candidate_parent = 0;

    for (size_t remote_i = 0; remote_i < sampled; ++remote_i)
    {
        if (deviate <= choose_prob[remote_i])
        {
            candidate_parent = remote_i;
            break;
        }
    }

    create_kid(
            MetaPop[patch_origin[candidate_parent]].locals[individual_ids[candidate_parent]],
            Kid);
    
    
    // keep track of new state variables
    // to update counters (see below)
    size_t n_adapted_new = n_adapted;

    // then mortality 
    assert(n_mortality_candidates > 0 && n_mortality_candidates <= Npp);

    // sample a random candidate for mortality
    size_t candidate_id = 
        mortality_candidate_ids[gsl_rng_uniform_int(r, n_mortality_candidates)];

    assert(candidate_id >= 0 && candidate_id < Npp);

    assert(
            (!mortality_maladapted &&
            MetaPop[random_patch_id].locals[candidate_id].z ==
            MetaPop[random_patch_id].state
            )
            ||
            (mortality_maladapted &&
            MetaPop[random_patch_id].locals[candidate_id].z !=
            MetaPop[random_patch_id].state
            )
    );

    if (!mortality_maladapted)
    {
        --n_adapted_new;
    }

    // mortality and replacement
    MetaPop[random_patch_id].locals[candidate_id] = Kid;

    if (MetaPop[random_patch_id].locals[candidate_id].z ==
            MetaPop[random_patch_id].state)
    {
        ++n_adapted_new;
    }

    // update statistics
    //
    // delete patch id in the corresponding stack
    // by replacing it with the last patch id in the stack
    // and then deleting the last element (by reducing the counter by 1)
    patch_ids[patch_envt][n_adapted][random_patch] =
        patch_ids[patch_envt][n_adapted][
            --Npatches[patch_envt][n_adapted] 
        ];

    // add patch id to the correct stack
    patch_ids[patch_envt][n_adapted_new][
        Npatches[patch_envt][n_adapted_new]++
    ] = random_patch_id;

    // we're done.
    // error checking
    assert(MetaPop[random_patch_id].state == patch_envt);

#ifndef NDEBUG
    n_adapted_check = 0;

    for (size_t breeder_i = 0; breeder_i < Npp; ++breeder_i)
    {
        n_adapted_check += MetaPop[random_patch_id].locals[breeder_i].z 
            == MetaPop[random_patch_id].state;
    }

    assert(n_adapted_check == n_adapted_new);
#endif
}

// main function body
int main(int argc, char **argv)
{
    init_arguments(argc, argv);
    init_pop();
    
    write_parameters();
    write_data_headers();

    // calculate probabilities of events 
    // there are a total number of 2 * (Npp + 1) * 3 events
    // i. because there are two environments
    // ii. each patch contains 0,1,2, .. Npp (=Npp + 1) 
    // combinations of adapted/maladapted indiviudals
    // iii. 3 events: 1. switching, 2. mortality adapted, 3. mortality maladapted
    size_t nprobs = 2 * (Npp + 1) * 3;

    double probs[nprobs]; // vector with probabilities

    bool done; // variable to get out of nested loop
    double prob; // actual event probability sampled from 
                // the cumulative prob distribution

    for (generation = 0; generation < NumGen; ++generation)
    {
        // generate cumulative prob. distribution of possible
        // events
        double sum_probs = 0;
        size_t prob_count = 0;

        // loop through the different patch types 
        for (size_t patch_envt = 0; patch_envt < 2; ++patch_envt)
        {
            // number of adapted individuals in a patch
            for (size_t n_adapted = 0; n_adapted <= Npp; ++n_adapted)
            {
                // probability that a patch switches state
                sum_probs += envt_switch[patch_envt] 
                    * Npatches[patch_envt][n_adapted];

                probs[prob_count++] = sum_probs; 

                // mortality adapted
                sum_probs += (double) n_adapted / Npp 
                    * C[patch_envt]
                    * Npatches[patch_envt][n_adapted];
                
                probs[prob_count++] = sum_probs; 

                // mortality maladapted individual
                sum_probs += (double) (Npp - n_adapted) / Npp 
                    * c[patch_envt] 
                    * Npatches[patch_envt][n_adapted];
                
                probs[prob_count++] = sum_probs; 

                assert(prob_count <= nprobs);
            }
        }

        // decide which event will happen
        //
        // draw value from cumul distribution
        prob = gsl_rng_uniform(r) *  sum_probs;

        prob_count = 0;

        // variable to break out of all the loops
        done = false;

        // locate the event
        // in the cumulative distribution
        for (size_t patch_envt = 0; patch_envt < 2; ++patch_envt)
        {
            if (done)
            {
                break;
            }

            for (size_t n_adapted = 0; n_adapted <= Npp; ++n_adapted)
            {
                if (prob <= probs[prob_count])
                {
                    switch_patch_state(patch_envt, n_adapted);
                    done = true;
                    break;
                }

                ++prob_count;

                // mortality adapted
                if (prob <= probs[prob_count])
                {
                    mortality(patch_envt, n_adapted, false);
                    done = true;
                    break;
                }

                ++prob_count;

                // mortality maladapted
                if (prob <= probs[prob_count])
                {
                    mortality(patch_envt, n_adapted, true);
                    done = true;
                    break;
                }
                
                ++prob_count;
            }
        } // end for patch_envt

        if (generation % output_interval == 0)
        {
            write_data();
        }
    }

    write_data();
}
