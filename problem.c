/**
 * Self-gravitating disc
 *
 * A self-gravitating disc is integrated using
 * the leap frog integrator. Collisions are not resolved.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "rebound.h"


void heartbeat(struct reb_simulation* const r);

void add_satellite(struct reb_simulation* const r, double inclination){
    // Satellite parameters
    double satellite_central_mass = 0.05; // Mass of satellite central body
    double satellite_disk_mass = 0.01;   // Mass of satellite disk
    int N_satellite = 500;              // Number of satellite disk particles

    // Physical size of satellite disk
    double a_min = 0.5;     // Minimum semi-major axis of satellite disk particles 
    double a_max = 1.8;     // Maximum semi-major axis of satellite disk particles

    // Initial orbital position
    double x0 = 8.0;
    double y0 = 0.0;
    double z0 = 0.0;

    // Total masses of primary and satellite
    double primary_mass = 1.0 + 0.2;
    double satellite_mass = satellite_central_mass + satellite_disk_mass;

    // Initial orbital distance and desired periapsis
    double r0 = x0;
    double rp = 2.5;

    // Total mass of the two galaxies
    double total_mass = primary_mass + satellite_mass;

    // Total speed for a parabolic orbit at r0
    double v_orbit = sqrt(
    2.0 * r->G * total_mass / r0
    );

    // Angular momentum for a parabolic orbit
    // with periapsis rp
    double h = sqrt(
        2.0 * r->G * total_mass * rp
    );

    // Tangential component of velocity
    double v_t = h / r0;

    // Radial component of velocity
    double v_r = sqrt(
        v_orbit*v_orbit - v_t*v_t
    );

    // Satellite initially moves inward toward the primary
    double vx_orbit = -v_r;

    // Incline the tangential velocity
    double vy_orbit = v_t * cos(inclination);
    double vz_orbit = v_t * sin(inclination);

    // Add satellite central mass
    struct reb_particle satellite_star = {0};

    satellite_star.m = satellite_central_mass;
    satellite_star.hash = 10001;

    satellite_star.x = x0;
    satellite_star.y = y0;
    satellite_star.z = z0;

    satellite_star.vx = vx_orbit;
    satellite_star.vy = vy_orbit;
    satellite_star.vz = vz_orbit;

    reb_simulation_add(r, satellite_star);

    // Add satellite disk particles
    for (int i = 0; i < N_satellite; i++){
        struct reb_particle pt = {0};
        pt.hash = 10002 + i;

        // Internal position within satellite disk
        double a = reb_random_powerlaw(
            r,
            a_min,
            a_max,
            -1.5
        );

        double phi = reb_random_uniform(
            r,
            0,
            2.0 * M_PI
        );

        pt.x = a * cos(phi) + x0;
        pt.y = a * sin(phi) + y0;
        pt.z = a * reb_random_normal(r, 0.001) + z0;

        // Satellite internal gravitational mass
        double mu = satellite_star.m
            + satellite_disk_mass *
            (
                pow(a, -3.0/2.0) - pow(a_min, -3.0/2.0)
            )
            /
            (
                pow(a_max, -3.0/2.0) - pow(a_min, -3.0/2.0)
            );

        // Internal circular velocity
        double vkep = sqrt(r->G * mu / a);

        // Internal disk velocity
        pt.vx = vkep * sin(phi) + vx_orbit;
        pt.vy = -vkep * cos(phi) + vy_orbit;
        pt.vz = vz_orbit;

        // Particle mass
        pt.m = satellite_disk_mass / (double)N_satellite;

        reb_simulation_add(r, pt);
    }
}

int main(int argc, char* argv[]){
    struct reb_simulation* const r = reb_simulation_create();
    
    r->rand_seed = 12345;   // Set the random seed for reproducibility

    // Start the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    // Setup constants
    r->integrator       = REB_INTEGRATOR_LEAPFROG; // Leapfrog integrator
    r->gravity          = REB_GRAVITY_TREE;        // Tree code gravity
    r->boundary         = REB_BOUNDARY_OPEN;       // Open boundary conditions
    r->opening_angle2   = 1.5;          // This constant determines the accuracy of the tree code gravity estimate.
    r->G                = 1;            // Gravitational constant
    r->softening        = 0.02;         // Gravitational softening length
    r->dt               = 1.5e-2;         // Timestep
    const double boxsize = 80.0;        // Size of the simulation box
    reb_simulation_configure_box(r,boxsize,1,1,1);  // Configure the simulation box for the tree code

    // Setup particles
    double disc_mass = 2e-1;    // Total disc mass
    int N = 10000;            // Number of particles
    // Initial conditions
    struct reb_particle star = {0};
    star.m         = 1;
    reb_simulation_add(r, star);
    for (int i=0;i<N;i++){
        struct reb_particle pt = {0};
        pt.hash = i + 1;
        double a_min = 1.02;
        double a_max = 4.25;
        double a = reb_random_powerlaw(r, a_min, a_max, -1.5);
        double phi     = reb_random_uniform(r, 0,2.*M_PI);
        pt.x         = a*cos(phi);
        pt.y         = a*sin(phi);
        pt.z         = a*reb_random_normal(r, 0.001);
        double mu = star.m + disc_mass *(pow(a,-3./2.)-pow(a_min,-3./2.)) /(pow(a_max,-3./2.)-pow(a_min,-3./2.));
        double vkep     = sqrt(r->G*mu/a);
        pt.vx         =  vkep * sin(phi);
        pt.vy         = -vkep * cos(phi);
        pt.vz         = 0;
        pt.m         = disc_mass/(double)N;
        reb_simulation_add(r, pt);
        
    }
    //Satellite Galaxy
    
    double inclination_deg = 0.0;       // Inclination of satellite orbit in degrees
    double inclination = inclination_deg * M_PI / 180.0;

    add_satellite(r, inclination);

    r->heartbeat = heartbeat;
    reb_simulation_integrate(r, INFINITY);
}

void heartbeat(struct reb_simulation* const r){
    if (reb_simulation_output_check(r, 1.0)){

        const int N_primary_disk = 10000;
        const double R_core = 2.5;
        const int N_core_min = 500;

        /*
         * ---------------------------------------------------------
         * 1. Find the primary-disk particles that are still
         *    present in the simulation.
         *
         *    Primary disk particles have hashes 1 through 10000.
         * ---------------------------------------------------------
         */

        int N_primary_present = 0;

        double total_mass = 0.0;
        double x_cm = 0.0;
        double y_cm = 0.0;
        double z_cm = 0.0;

        for (int i = 0; i < r->N; i++){

            struct reb_particle p = r->particles[i];

            if (p.hash >= 1 && p.hash <= N_primary_disk){

                N_primary_present++;

                total_mass += p.m;

                x_cm += p.m * p.x;
                y_cm += p.m * p.y;
                z_cm += p.m * p.z;
            }
        }

        /*
         * If somehow no primary particles remain, stop the
         * diagnostic calculation.
         */

        if (N_primary_present == 0){
            printf("%f 0 0 0 0 0 0 0 0 0 0 0 0\n", r->t);
            return;
        }

        x_cm /= total_mass;
        y_cm /= total_mass;
        z_cm /= total_mass;


        /*
         * ---------------------------------------------------------
         * 2. Diagnostics for the entire primary disk
         * ---------------------------------------------------------
         */

        double R2_sum_all = 0.0;
        double z2_sum_all = 0.0;
        double Lz_all = 0.0;

        for (int i = 0; i < r->N; i++){

            struct reb_particle p = r->particles[i];

            if (p.hash >= 1 && p.hash <= N_primary_disk){

                double dx = p.x - x_cm;
                double dy = p.y - y_cm;
                double dz = p.z - z_cm;

                double R2 = dx*dx + dy*dy;

                R2_sum_all += R2;
                z2_sum_all += dz*dz;

                Lz_all += p.m * (
                    dx*p.vy - dy*p.vx
                );
            }
        }

        double R_rms_all = sqrt(
            R2_sum_all / N_primary_present
        );

        double z_rms_all = sqrt(
            z2_sum_all / N_primary_present
        );


        /*
         * ---------------------------------------------------------
         * 3. Find the primary particles inside the core
         *
         *    Core definition:
         *
         *                 R < 2.5
         *
         *    R is measured relative to the primary disk COM.
         * ---------------------------------------------------------
         */

        int N_core = 0;

        double core_mass = 0.0;

        double x_core_cm = 0.0;
        double y_core_cm = 0.0;
        double z_core_cm = 0.0;

        for (int i = 0; i < r->N; i++){

            struct reb_particle p = r->particles[i];

            if (p.hash >= 1 && p.hash <= N_primary_disk){

                double dx = p.x - x_cm;
                double dy = p.y - y_cm;

                double R2 = dx*dx + dy*dy;

                if (R2 < R_core*R_core){

                    N_core++;

                    core_mass += p.m;

                    x_core_cm += p.m*p.x;
                    y_core_cm += p.m*p.y;
                    z_core_cm += p.m*p.z;
                }
            }
        }


        /*
         * ---------------------------------------------------------
         * 4. Core centre of mass
         * ---------------------------------------------------------
         */

        if (N_core > 0){

            x_core_cm /= core_mass;
            y_core_cm /= core_mass;
            z_core_cm /= core_mass;
        }


        /*
         * ---------------------------------------------------------
         * 5. Core RMS quantities
         * ---------------------------------------------------------
         */

        double R2_sum_core = 0.0;
        double z2_sum_core = 0.0;

        if (N_core >= N_core_min){

            for (int i = 0; i < r->N; i++){

                struct reb_particle p = r->particles[i];

                if (p.hash >= 1 && p.hash <= N_primary_disk){

                    double dx_global = p.x - x_cm;
                    double dy_global = p.y - y_cm;

                    double R2_global =
                        dx_global*dx_global +
                        dy_global*dy_global;

                    if (R2_global < R_core*R_core){

                        double dx_core =
                            p.x - x_core_cm;

                        double dy_core =
                            p.y - y_core_cm;

                        double dz_core =
                            p.z - z_core_cm;

                        double R2_core =
                            dx_core*dx_core +
                            dy_core*dy_core;

                        R2_sum_core += R2_core;
                        z2_sum_core += dz_core*dz_core;
                    }
                }
            }
        }


        /*
         * ---------------------------------------------------------
         * 6. Core RMS values
         *
         *    If fewer than 500 primary particles are in the
         *    core, mark the measurement as invalid with -1.
         * ---------------------------------------------------------
         */

        double R_rms_core = -1.0;
        double z_rms_core = -1.0;

        if (N_core >= N_core_min){

            R_rms_core = sqrt(
                R2_sum_core / N_core
            );

            z_rms_core = sqrt(
                z2_sum_core / N_core
            );
        }


        /*
         * ---------------------------------------------------------
         * 7. Fraction of original primary disk inside core
         * ---------------------------------------------------------
         */

        double f_core =
            (double)N_core /
            (double)N_primary_disk;


        /*
         * ---------------------------------------------------------
         * 8. Fraction of original primary disk still present
         * ---------------------------------------------------------
         */

        double f_primary_present =
            (double)N_primary_present /
            (double)N_primary_disk;


        /*
         * ---------------------------------------------------------
         * 9. Total energy
         * ---------------------------------------------------------
         */

        double E = reb_simulation_energy(r);


        /*
         * ---------------------------------------------------------
         * 10. Output
         *
         * Columns:
         *
         *  1. t
         *  2. N_primary_present
         *  3. f_primary_present
         *  4. N_core
         *  5. f_core
         *  6. R_rms_all
         *  7. z_rms_all
         *  8. R_rms_core
         *  9. z_rms_core
         * 10. z_core_cm
         * 11. E
         * 12. Lz_all
         * 13. x_cm
         * 14. z_cm
         * ---------------------------------------------------------
         */

        printf(
            "%f %d %e %d %e %e %e %e %e %e %e %e %e %e\n",
            r->t,
            N_primary_present,
            f_primary_present,
            N_core,
            f_core,
            R_rms_all,
            z_rms_all,
            R_rms_core,
            z_rms_core,
            z_core_cm,
            E,
            Lz_all,
            x_cm,
            z_cm
        );
    }
}