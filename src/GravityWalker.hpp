/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2015 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * Shadowfax is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Shadowfax is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with Shadowfax. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file GravityWalker.hpp
 *
 * @brief Treewalkers used for gravity tree walks: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEAD_GRAVITYWALKER
#define HEAD_GRAVITYWALKER

#include "MPIMethods.hpp"               // for MyMPI_Pack, MyMPI_Unpack
#include "utilities/GasParticle.hpp"    // for GasParticle
#include "utilities/Particle.hpp"       // for Particle
#include "utilities/ParticleTypes.hpp"  // for ParticleType, etc
#include "utilities/TreeWalker.hpp"     // for PeriodicTreeWalker
#include <stddef.h>                     // for NULL

class EwaldTable;
class Leaf;
class PseudoNode;
class TreeNode;

/**
 * @brief Interface for gravity kernels
 *
 * Defines a functor that gives the kernel values, and functions that give its
 * primitive and derivative.
 */
class GravityKernel {
  public:
    virtual ~GravityKernel() {}

    /**
     * @brief Evaluate the kernel at the given relative radius and with the
     * given inverse kernel volume
     *
     * @param u Relative radius (distance in units of softening length)
     * @param hsoftinv3 Inverse kernel volume (inverse softening length cubed)
     * @return Value of the kernel at the given relative radius
     */
    virtual double operator()(double u, double hsoftinv3) = 0;

    /**
     * @brief Get the value of the integral of the kernel at the given inverse
     * radius
     *
     * @param u Relative radius (distance in units of softening length)
     * @param hsoftinv Inverse softening length
     * @return Value of the integral of the kernel at the given relative radius
     */
    virtual double primitive(double u, double hsoftinv) = 0;

    /**
     * @brief Get the value of the derivative of the kernel at the given inverse
     * radius
     *
     * @param u Relative radius (distance in units of softening length)
     * @param hsoftinv Inverse softening length
     * @return Value of the derivative of the kernel at the given relative
     * radius
     */
    virtual double derivative(double u, double hsoftinv) = 0;
};

/**
 * @brief Softening kernel used for GasParticles.
 */
class GasGravityKernel : public GravityKernel {
  public:
    /**
     * @brief Top-hat kernel from AREPO paper
     *
     * @param u Distance divided by softening length
     * @param hsoftinv3 Inverse softening length cubed
     * @return Value of hte kernel at the given distance
     */
    virtual double operator()(double u, double hsoftinv3) {
        return hsoftinv3;
    }

    /**
     * @brief Get the value of the integral of the kernel at the given inverse
     * radius
     *
     * @param u Relative radius (distance in units of softening length)
     * @param hsoftinv Inverse softening length
     * @return Value of the integral of the kernel at the given relative radius
     */
    virtual double primitive(double u, double hsoftinv) {
        return -hsoftinv * 0.5 * (3. - u * u);
    }

    /**
     * @brief Get the value of the derivative of the kernel at the given inverse
     * radius
     *
     * @param u Relative radius (distance in units of softening length)
     * @param hsoftinv Inverse softening length
     * @return Value of the derivative of the kernel at the given relative
     * radius
     */
    virtual double derivative(double u, double hsoftinv) {
        double hsoftinv2 = hsoftinv * hsoftinv;
        return 1.5 * (hsoftinv2 - hsoftinv2 * hsoftinv2);
    }
};

/**
 * @brief GravityKernel used for DM particles
 */
class DMGravityKernel : public GravityKernel {
  public:
    /**
     * @brief The Gadget-2 gravity kernel
     *
     * @param u Distance divided by softening length
     * @param hsoftinv3 Inverse softening length cubed
     * @return Value of the kernel at the given distance
     */
    virtual double operator()(double u, double hsoftinv3) {
        double u2 = u * u;
        if(u < 0.5) {
            return hsoftinv3 * (10.666666666667 + u2 * (32.0 * u - 38.4));
        } else {
            double u3 = u * u2;
            return hsoftinv3 * (21.333333333333 - 48.0 * u + 38.4 * u2 -
                                10.666666666667 * u3 - 0.066666666667 / (u3));
        }
    }

    /**
     * @brief Get the value of the integral of the kernel at the given inverse
     * radius
     *
     * @param u Relative radius (distance in units of softening length)
     * @param hsoftinv Inverse softening length
     * @return Value of the integral of the kernel at the given relative radius
     */
    virtual double primitive(double u, double hsoftinv) {
        double u2 = u * u;
        double kernel;
        if(u < 0.5) {
            kernel = -2.8 + u2 * (5.333333333333 + u2 * (6.4 * u - 9.6));
        } else {
            kernel = -3.2 + 0.066666666667 / u +
                     u2 * (10.666666666667 +
                           u * (-16.0 + u * (9.6 - 2.133333333333 * u)));
        }
        return hsoftinv * kernel;
    }

    /**
     * @brief Get the value of the derivative of the kernel at the given inverse
     * radius
     *
     * @param u Relative radius (distance in units of softening length)
     * @param hsoftinv Inverse softening length
     * @return Value of the derivative of the kernel at the given relative
     * radius
     */
    virtual double derivative(double u, double hsoftinv) {
        return 0.;
    }
};

/**
 * @brief TreeWalker to calculate the gravitational acceleration using a
 * relative opening criterion as described in Springel 2005
 */
class GravityWalker : public PeriodicTreeWalker {
  private:
    /*! @brief Particle for which the acceleration is calculated */
    Particle* _p;

    /*! @brief Position of the Particle for which the acceleration is
     *  calculated */
    Vec _position;

    /*! @brief Current value of the gravitational acceleration */
    Vec _a;

    /*! @brief Size of the previous gravitational acceleration, used for the
     *  relative opening criterion */
    double _olda;

    /*! @brief Softening length of the Particle */
    double _hsoft;

    /*! @brief Inverse softening length of the Particle */
    double _hsoftinv;

    /*! @brief Inverse softening length cubed */
    double _hsoftinv3;

    /*! @brief Flag specifying whether the treewalk is local or is being done
     *  for a Particle residing on another MPI process */
    bool _local;

    /*! @brief Timer for the computational cost of this treewalk */
    Timer _comp_cost;

    /*! @brief GravityKernel used to calculate softened accelerations */
    GravityKernel* _kernel;

    /*! @brief \f$\eta\f$-parameter for variable softening lengths */
    double _eta;

    /*! @brief Type of the Particle (PARTTYPE_GAS/PARTTYPE_DM) */
    ParticleType _type;

  public:
    /**
     * @brief Auxiliary class used to communicate particle information between
     * MPI processes during the gravity treewalk.
     */
    class Export {
      private:
        /*! @brief Particle for which the treewalk is performed */
        Particle* _p;

        /*! @brief Size of the previous gravitational acceleration, used for the
         *  relative opening criterion (will be exported) */
        double _olda;

        /*! @brief Position of the Particle (will be exported to the other
         *  process) */
        Vec _pos;

        /*! @brief Softening length of the Particle (will be exported to the
         *  other process) */
        double _hsoft;

        /*! @brief Type of the Particle (PARTTYPE_GAS/PARTTYPE_DM) (will be
         *  exported) */
        ParticleType _type;

      public:
        /**
         * @brief Constructor
         *
         * @param p Particle for which the treewalk is performed
         * @param olda Old acceleration of the Particle
         */
        Export(Particle* p, double olda) {
            _p = p;
            _olda = olda;
            _pos = p->get_position();
            _hsoft = p->get_hsoft();
            _type = p->type();
        }

        /**
         * @brief Add the relevant particle data to the given communication
         * buffer
         *
         * @param buffer Buffer to write to
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            MyMPI_Pack(&_pos[0], ndim_, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_olda, 1, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_hsoft, 1, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_type, 1, MPI_INT, buffer, bufsize, position);
        }

        /**
         * @brief Read the response of the external communication from the given
         * communication buffer
         *
         * @param buffer Buffer to read from
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void unpack_data(void* buffer, int bufsize, int* position) {
            Vec a;
            float comp_cost;
            MyMPI_Unpack(buffer, bufsize, position, &a[0], ndim_, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &comp_cost, 1, MPI_FLOAT);
            _p->set_gravitational_acceleration(
                    _p->get_gravitational_acceleration() + a);
            _p->add_comp_cost(comp_cost);
            if(_type == PARTTYPE_GAS) {
                double eta;
                MyMPI_Unpack(buffer, bufsize, position, &eta, 1, MPI_DOUBLE);
                GasParticle* gas = (GasParticle*)_p;
                gas->set_eta(gas->get_eta() + eta);
            }
        }
    };

    /**
     * @brief Auxiliary class used to communicate particle information between
     * MPI processes during the gravity treewalk.
     */
    class Import {
      private:
        /*! @brief Position for which the treewalk is performed */
        Vec _pos;

        /*! @brief Old acceleration, used for the relative opening criterion */
        double _olda;

        /*! @brief Local result of the treewalk (is exported) */
        Vec _a;

        /*! @brief Computational cost of the local treewalk (is exported) */
        float _comp_cost;

        /*! @brief Softening length of the Particle for which the treewalk is
         *  performed */
        double _hsoft;

        /*! @brief Type of the Particle for which the treewalk is performed
         *  (PARTTYPE_GAS/PARTTYPE_DM) */
        ParticleType _type;

        /*! @brief \f$\eta\f$-parameter for variable softening lengths (is
            exported) */
        double _eta;

      public:
        /**
         * @brief Constructor. Initialize the import based on the given
         * communication stream
         *
         * @param buffer Buffer to read from
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        Import(void* buffer, int bufsize, int* position) {
            MyMPI_Unpack(buffer, bufsize, position, &_pos[0], ndim_,
                         MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_olda, 1, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_hsoft, 1, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_type, 1, MPI_INT);
        }

        /**
         * @brief Get the position for which the treewalk is performed
         *
         * @return The position of the external Particle
         */
        Vec get_position() {
            return _pos;
        }

        /**
         * @brief Get the old acceleration of the external Particle
         *
         * @return Old acceleration of the external Particle
         */
        double get_olda() {
            return _olda;
        }

        /**
         * @brief Get the softening length of the external Particle
         *
         * @return Softening length of the external Particle
         */
        double get_hsoft() {
            return _hsoft;
        }

        /**
         * @brief Get the type of the external Particle
         *
         * @return Type of the external Particle (PARTTYPE_GAS/PARTTYPE_DM)
         */
        ParticleType get_type() {
            return _type;
        }

        /**
         * @brief set_a Set the acceleration
         *
         * @param a Value of the acceleration
         */
        void set_a(Vec a) {
            _a = a;
        }

        /**
         * @brief Set the computational cost
         * @param comp_cost Value of the computational cost
         */
        void set_comp_cost(float comp_cost) {
            _comp_cost = comp_cost;
        }

        /**
         * @brief Set the \f$\eta\f$-factor for variable softening lengths
         *
         * @param eta Value of the \f$\eta\f$-factor
         */
        void set_eta(double eta) {
            _eta = eta;
        }

        /**
         * @brief Add the relevant data to the given communication buffer to
         * send them back to the original process
         *
         * @param buffer Communication buffer
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            MyMPI_Pack(&_a[0], ndim_, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_comp_cost, 1, MPI_FLOAT, buffer, bufsize, position);
            if(_type == PARTTYPE_GAS) {
                MyMPI_Pack(&_eta, 1, MPI_DOUBLE, buffer, bufsize, position);
            }
        }
    };

    GravityWalker(Particle* p, bool local = true);
    GravityWalker(Import& import);

    ~GravityWalker();

    Vec get_acceleration();

    bool splitnode(TreeNode* node);
    void nodeaction(TreeNode* node);
    void leafaction(Leaf* leaf);
    void pseudonodeaction(PseudoNode* pseudonode);
    bool export_to_pseudonode(PseudoNode* pseudonode);

    void set_position(Vec position);
    Vec get_position();

    bool periodicsplitnode(TreeNode* node, EwaldTable& ewald_table);
    void periodicpseudonodeaction(PseudoNode* pseudonode,
                                  EwaldTable& ewald_table);
    void periodicleafaction(Leaf* leaf, EwaldTable& ewald_table);

    void after_walk();

    void after_walk(Import& import);

    Export get_export();
};

/**
 * @brief TreeWalker used to calculate the Ewald correction to the gravitational
 * force
 */
class EwaldGravityWalker : public PeriodicTreeWalker {
  private:
    /*! @brief Particle for which the Ewald correction force is calculated */
    Particle* _p;

    /*! @brief Position of the Particle for which the Ewald correction force is
     *  calculated */
    Vec _position;

    /*! @brief Current value of the Ewald correction force */
    Vec _acorr;

    /*! @brief Reference to the EwaldTable used to calculate corrections */
    EwaldTable& _ewald_table;

    /*! @brief Old gravitational acceleration */
    double _olda;

  public:
    /**
     * @brief Auxiliary class used to communicate particle information to other
     * MPI processes during the Ewald correction treewalk
     */
    class Export {
      private:
        /*! @brief Particle for which the Ewald correction is calculated */
        Particle* _p;

        /*! @brief Old acceleration of the Particle (is exported) */
        double _olda;

        /*! @brief Position of the Particle (is exported) */
        Vec _pos;

        /*! @brief Softening length of the Particle (is exported) */
        double _hsoft;

        /*! @brief Type of the Particle (PARTTYPE_GAS/PARTTYPE_DM) (is
         *  exported) */
        ParticleType _type;

      public:
        /**
         * @brief Constructor
         *
         * @param p Particle for which the Ewald correction is calculated
         * @param olda Old acceleration of the Particle
         */
        Export(Particle* p, double olda) {
            _p = p;
            _olda = olda;
            _pos = p->get_position();
            _hsoft = p->get_hsoft();
            _type = p->type();
        }

        /**
         * @brief Add relevant data to the given communication buffer for export
         * to another MPI-process
         *
         * @param buffer Communication buffer
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            MyMPI_Pack(&_pos[0], ndim_, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_olda, 1, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_hsoft, 1, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_type, 1, MPI_INT, buffer, bufsize, position);
        }

        /**
         * @brief Read the response from the external treewalk from the given
         * communication buffer and finalize the external treewalk
         *
         * @param buffer Communication buffer
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void unpack_data(void* buffer, int bufsize, int* position) {
            Vec acorr;
            float comp_cost;
            MyMPI_Unpack(buffer, bufsize, position, &acorr[0], ndim_,
                         MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &comp_cost, 1, MPI_FLOAT);
            _p->set_gravitational_acceleration(
                    _p->get_gravitational_acceleration() + acorr);
            _p->add_comp_cost(comp_cost);
        }
    };

    /**
     * @brief Auxiliary class used to communicate particle information to other
     * MPI processes during the Ewald correction treewalk
     */
    class Import {
      private:
        /*! @brief Position for which the Ewald correction is calculated */
        Vec _pos;

        /*! @brief Old acceleration, used for the relative opening criterion */
        double _olda;

        /*! @brief Result of the local treewalk (is exported) */
        Vec _acorr;

        /*! @brief Computational cost of the local treewalk (is exported) */
        float _comp_cost;

        /*! @brief Softening length of the external Particle */
        double _hsoft;

        /*! @brief Type of the external Particle (PARTTYPE_GAS/PARTTYPE_DM) */
        ParticleType _type;

      public:
        /**
         * @brief Constructor. Initialize the import from the given
         * communication buffer
         *
         * @param buffer Communication buffer
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        Import(void* buffer, int bufsize, int* position) {
            MyMPI_Unpack(buffer, bufsize, position, &_pos[0], ndim_,
                         MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_olda, 1, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_hsoft, 1, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_type, 1, MPI_INT);
        }

        /**
         * @brief Get the position of the external Particle
         *
         * @return Position of the external Particle
         */
        Vec get_position() {
            return _pos;
        }

        /**
         * @brief Get the old acceleration of the external Particle
         *
         * @return Old acceleration of the Particle
         */
        double get_olda() {
            return _olda;
        }

        /**
         * @brief Get the softening length of the external Particle
         *
         * @return Softening length of the external Particle
         */
        double get_hsoft() {
            return _hsoft;
        }

        /**
         * @brief Get the type of the external Particle
         *
         * @return Type of the external Particle (PARTTYPE_GAS/PARTTYPE_DM)
         */
        ParticleType get_type() {
            return _type;
        }

        /**
         * @brief Set the Ewald correction force
         *
         * @param acorr Value of the Ewald correction
         */
        void set_acorr(Vec acorr) {
            _acorr = acorr;
        }

        /**
         * @brief Set the computational cost
         *
         * @param comp_cost Value of the computational cost
         */
        void set_comp_cost(float comp_cost) {
            _comp_cost = comp_cost;
        }

        /**
         * @brief Add relevant data to the given communication stream to send
         * them back to the original process
         *
         * @param buffer Communication buffer
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            MyMPI_Pack(&_acorr[0], ndim_, MPI_DOUBLE, buffer, bufsize,
                       position);
            MyMPI_Pack(&_comp_cost, 1, MPI_FLOAT, buffer, bufsize, position);
        }
    };

    EwaldGravityWalker(Particle* p, EwaldTable& ewald_table, bool local = true);
    EwaldGravityWalker(Import& import, EwaldTable& ewald_table);

    ~EwaldGravityWalker();

    Vec get_acceleration();

    bool splitnode(TreeNode* node);
    void nodeaction(TreeNode* node);
    void leafaction(Leaf* leaf);
    void pseudonodeaction(PseudoNode* pseudonode);
    bool export_to_pseudonode(PseudoNode* pseudonode);

    void set_position(Vec position);
    Vec get_position();

    bool periodicsplitnode(TreeNode* node, EwaldTable& ewald_table);
    void periodicpseudonodeaction(PseudoNode* pseudonode,
                                  EwaldTable& ewald_table);
    void periodicleafaction(Leaf* leaf, EwaldTable& ewald_table);

    void after_walk();

    void after_walk(Import& import);

    Export get_export();
};

/**
 * @brief TreeWalker used to calculate the gravitational potential for a given
 * Particle
 */
class PotentialWalker : public PeriodicTreeWalker {
  private:
    /*! @brief Particle for which the potential is calculated */
    Particle* _p;

    /*! @brief Position of the Particle for which the potential is calculated */
    Vec _position;

    /*! @brief Current value of the potential */
    double _epot;

    /*! @brief Old gravitational acceleration, used for the relative opening
     *  criterion */
    double _olda;

    /*! @brief Softening length of the Particle for which the potential is
     *  calculated */
    double _hsoft;

    /*! @brief Inverse softening length of the Particle for which the potential
     *  is calculated */
    double _hsoftinv;

    /*! @brief Inverse softening length cubed */
    double _hsoftinv3;

    /*! @brief Flag indicating if the treewalk is performed for a local
     *  Particle */
    bool _local;

    /*! @brief Timer for the computational cost of the treewalk */
    Timer _comp_cost;

    /*! @brief GravityKernel used for close interaction softening */
    GravityKernel* _kernel;

  public:
    /**
     * @brief Auxiliary class used to communicate particle information to other
     * MPI processes during the gravity potential treewalk
     */
    class Export {
      private:
        /*! @brief Particle for which the potential is calculated */
        Particle* _p;

        /*! @brief Old acceleration of the Particle (is exported) */
        double _olda;

        /*! @brief Position of the Particle (is exported) */
        Vec _pos;

        /*! @brief Softening length of the Particle (is exported) */
        double _hsoft;

        /*! @brief Type of the Particle (PARTTYPE_GAS/PARTTYPE_DM) (is
         *  exported) */
        ParticleType _type;

      public:
        /**
         * @brief Constructor
         *
         * @param p Particle for which the potential is calculated
         * @param olda Old acceleration of the Particle
         */
        Export(Particle* p, double olda) {
            _p = p;
            _olda = olda;
            _pos = p->get_position();
            _hsoft = p->get_hsoft();
            _type = p->type();
        }

        /**
         * @brief Add relevant data to the given communication buffer for export
         * to another MPI-process
         *
         * @param buffer Communication buffer
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            MyMPI_Pack(&_pos[0], ndim_, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_olda, 1, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_hsoft, 1, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_type, 1, MPI_INT, buffer, bufsize, position);
        }

        /**
         * @brief Read the response from the external treewalk from the given
         * communication buffer and finalize the external treewalk
         *
         * @param buffer Communication buffer
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void unpack_data(void* buffer, int bufsize, int* position) {
            double epot;
            float comp_cost;
            MyMPI_Unpack(buffer, bufsize, position, &epot, 1, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &comp_cost, 1, MPI_FLOAT);
            _p->set_gravitational_potential(_p->get_gravitational_potential() +
                                            epot);
            _p->add_comp_cost(comp_cost);
        }
    };

    /**
     * @brief Auxiliary class used to communicate particle information to other
     * MPI processes during the gravity potential treewalk
     */
    class Import {
      private:
        /*! @brief Position for which the potential is calculated */
        Vec _pos;

        /*! @brief Old acceleration, used for the relative opening criterion */
        double _olda;

        /*! @brief Result of the local treewalk (is exported) */
        double _epot;

        /*! @brief Computational cost of the local treewalk (is exported) */
        float _comp_cost;

        /*! @brief Softening length of the external Particle */
        double _hsoft;

        /*! @brief Type of the external Particle (PARTTYPE_GAS/PARTTYPE_DM) */
        ParticleType _type;

      public:
        /**
         * @brief Constructor. Initialize the import from the given
         * communication buffer
         *
         * @param buffer Communication buffer
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        Import(void* buffer, int bufsize, int* position) {
            MyMPI_Unpack(buffer, bufsize, position, &_pos[0], ndim_,
                         MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_olda, 1, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_hsoft, 1, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_type, 1, MPI_INT);
        }

        /**
         * @brief Get the position of the external Particle
         *
         * @return Position of the external Particle
         */
        Vec get_position() {
            return _pos;
        }

        /**
         * @brief Get the old acceleration of the external Particle
         *
         * @return Old acceleration of the Particle
         */
        double get_olda() {
            return _olda;
        }

        /**
         * @brief Get the softening length of the external Particle
         *
         * @return Softening length of the external Particle
         */
        double get_hsoft() {
            return _hsoft;
        }

        /**
         * @brief Get the type of the external Particle
         *
         * @return Type of the external Particle (PARTTYPE_GAS/PARTTYPE_DM)
         */
        ParticleType get_type() {
            return _type;
        }

        /**
         * @brief Set the potential
         *
         * @param epot Value of the gravitational potential
         */
        void set_epot(double epot) {
            _epot = epot;
        }

        /**
         * @brief Set the computational cost
         *
         * @param comp_cost Value of the computational cost
         */
        void set_comp_cost(float comp_cost) {
            _comp_cost = comp_cost;
        }

        /**
         * @brief Add relevant data to the given communication stream to send
         * them back to the original process
         *
         * @param buffer Communication buffer
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            MyMPI_Pack(&_epot, 1, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_comp_cost, 1, MPI_FLOAT, buffer, bufsize, position);
        }
    };

    PotentialWalker(Particle* p, bool local = true);
    PotentialWalker(Import& import);

    ~PotentialWalker();

    double get_epot();

    bool splitnode(TreeNode* node);
    void nodeaction(TreeNode* node);
    void leafaction(Leaf* leaf);
    void pseudonodeaction(PseudoNode* pseudonode);
    bool export_to_pseudonode(PseudoNode* pseudonode);

    void set_position(Vec position);
    Vec get_position();

    bool periodicsplitnode(TreeNode* node, EwaldTable& ewald_table);
    void periodicpseudonodeaction(PseudoNode* pseudonode,
                                  EwaldTable& ewald_table);
    void periodicleafaction(Leaf* leaf, EwaldTable& ewald_table);

    void after_walk();

    void after_walk(Import& import);

    Export get_export();
};

/**
 * @brief TreeWalker to perform a gravity force treewalk with a Barnes-Hut
 * opening criterion
 */
class BHGravityWalker : public PeriodicTreeWalker {
  private:
    /*! @brief Particle for which the gravitational acceleration is
     *  calculated */
    Particle* _p;

    /*! @brief Position of the Particle for which the gravitational acceleration
     *  is calculated */
    Vec _position;

    /*! @brief Current value of the gravitational acceleration */
    Vec _a;

    /*! @brief Softening lenght of the Particle */
    double _hsoft;

    /*! @brief Inverse softening length of the Particle */
    double _hsoftinv;

    /*! @brief Inverse softening length of the Particle cubed */
    double _hsoftinv3;

    /*! @brief Flag indicating whether the treewalk is performed for a local
     *  Particle */
    bool _local;

    /*! @brief Timer for the computational cost of the treewalk */
    Timer _comp_cost;

    /*! @brief GravityKernel used for close interaction softening */
    GravityKernel* _kernel;

  public:
    /**
     * @brief Auxiliary class used to communicate particle information to other
     * MPI processes during the BH gravity treewalk
     */
    class Export {
      private:
        /*! @brief Particle for which the gravitational acceleration is
         *  calculated */
        Particle* _p;

        /*! @brief Position of the Particle (is exported) */
        Vec _pos;

        /*! @brief Softening length of the Particle (is exported) */
        double _hsoft;

        /*! @brief Type of the Particle (PARTTYPE_GAS/PARTTYPE_DM) (is
         *  exported) */
        ParticleType _type;

      public:
        /**
         * @brief Constructor
         *
         * @param p Particle for which the gravitational acceleration is
         * calculated
         */
        Export(Particle* p) {
            _p = p;
            _pos = p->get_position();
            _hsoft = p->get_hsoft();
        }

        /**
         * @brief Add relevant data to the given communication buffer for export
         * to another MPI-process
         *
         * @param buffer Communication buffer
         * @param bufsize Buffer size
         * @param position Current position of the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            MyMPI_Pack(&_pos[0], ndim_, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_hsoft, 1, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_type, 1, MPI_INT, buffer, bufsize, position);
        }

        /**
         * @brief Read the response from an external process from the given
         * communication buffer
         *
         * @param buffer Communication buffer
         * @param bufsize Buffer size
         * @param position Current position of the buffer (is updated)
         */
        void unpack_data(void* buffer, int bufsize, int* position) {
            Vec a;
            float comp_cost;
            MyMPI_Unpack(buffer, bufsize, position, &a[0], ndim_, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &comp_cost, 1, MPI_FLOAT);
            _p->set_gravitational_acceleration(
                    _p->get_gravitational_acceleration() + a);
            _p->add_comp_cost(comp_cost);
        }
    };

    /**
     * @brief Auxiliary class used to communicate particle information during
     * the BH gravity treewalk
     */
    class Import {
      private:
        /*! @brief Position of the external Particle */
        Vec _pos;

        /*! @brief Result of the local treewalk (is exported) */
        Vec _a;

        /*! @brief Computational cost of the local treewalk (is exported) */
        float _comp_cost;

        /*! @brief Softening length of the external Particle */
        double _hsoft;

        /*! @brief Type of the external Particle (PARTTYPE_GAS/PARTTYPE_DM) */
        ParticleType _type;

      public:
        /**
         * @brief Constructor. Initialize the import based on the given
         * communication buffer
         *
         * @param buffer Communication buffer
         * @param bufsize Buffer size
         * @param position Current position of the buffer (is updated)
         */
        Import(void* buffer, int bufsize, int* position) {
            MyMPI_Unpack(buffer, bufsize, position, &_pos[0], ndim_,
                         MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_hsoft, 1, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_type, 1, MPI_INT);
        }

        /**
         * @brief Get the position of the external Particle
         *
         * @return Position of the external Particle
         */
        Vec get_position() {
            return _pos;
        }

        /**
         * @brief Get the softening length of the external Particle
         *
         * @return Softening length of the external Particle
         */
        double get_hsoft() {
            return _hsoft;
        }

        /**
         * @brief Get the type of the external Particle
         *
         * @return Type of the external Particle (PARTTYPE_GAS/PARTTYPE_DM)
         */
        ParticleType get_type() {
            return _type;
        }

        /**
         * @brief Set the gravitational acceleration
         *
         * @param a Value of the gravitational acceleration
         */
        void set_a(Vec a) {
            _a = a;
        }

        /**
         * @brief Set the computational cost
         *
         * @param comp_cost Value of the computational cost
         */
        void set_comp_cost(float comp_cost) {
            _comp_cost = comp_cost;
        }

        /**
         * @brief Add relevant data to the given communication buffer to send
         * them back to the original process
         *
         * @param buffer Communication buffer
         * @param bufsize Buffer size
         * @param position Current position of the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            MyMPI_Pack(&_a[0], ndim_, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_comp_cost, 1, MPI_FLOAT, buffer, bufsize, position);
        }
    };

    BHGravityWalker(Particle* p, bool local = true);
    BHGravityWalker(Import& import);

    ~BHGravityWalker();

    Vec get_acceleration();

    bool splitnode(TreeNode* node);
    void nodeaction(TreeNode* node);
    void leafaction(Leaf* leaf);
    void pseudonodeaction(PseudoNode* pseudonode);
    bool export_to_pseudonode(PseudoNode* pseudonode);

    void set_position(Vec position);
    Vec get_position();

    bool periodicsplitnode(TreeNode* node, EwaldTable& ewald_table);
    void periodicpseudonodeaction(PseudoNode* pseudonode,
                                  EwaldTable& ewald_table);
    void periodicleafaction(Leaf* leaf, EwaldTable& ewald_table);

    void after_walk();

    void after_walk(Import& import);

    Export get_export();
};

/**
 * @brief TreeWalker used to calculate the Ewald correction to the gravitational
 * force using a Barnes-Hut opening criterion
 */
class BHEwaldGravityWalker : public PeriodicTreeWalker {
  private:
    /*! @brief Particle for which the Ewald correction force is calculated */
    Particle* _p;

    /*! @brief Position of the Particle for which the Ewald correction force is
     *  calculated */
    Vec _position;

    /*! @brief Current value of the Ewald correction force */
    Vec _acorr;

    /*! @brief Reference to the EwaldTable used to calculate corrections */
    EwaldTable& _ewald_table;

  public:
    /**
     * @brief Auxiliary class used to communicate particle information to other
     * MPI processes during the Ewald correction treewalk
     */
    class Export {
      private:
        /*! @brief Particle for which the Ewald correction is calculated */
        Particle* _p;

        /*! @brief Old acceleration of the Particle (is exported) */
        double _olda;

        /*! @brief Position of the Particle (is exported) */
        Vec _pos;

        /*! @brief Softening length of the Particle (is exported) */
        double _hsoft;

        /*! @brief Type of the Particle (PARTTYPE_GAS/PARTTYPE_DM) (is
         *  exported) */
        ParticleType _type;

      public:
        /**
         * @brief Constructor
         *
         * @param p Particle for which the Ewald correction is calculated
         * @param olda Old acceleration of the Particle
         */
        Export(Particle* p, double olda) {
            _p = p;
            _olda = olda;
            _pos = p->get_position();
            _hsoft = p->get_hsoft();
            _type = p->type();
        }

        /**
         * @brief Add relevant data to the given communication buffer for export
         * to another MPI-process
         *
         * @param buffer Communication buffer
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            MyMPI_Pack(&_pos[0], ndim_, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_olda, 1, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_hsoft, 1, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_type, 1, MPI_INT, buffer, bufsize, position);
        }

        /**
         * @brief Read the response from the external treewalk from the given
         * communication buffer and finalize the external treewalk
         *
         * @param buffer Communication buffer
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void unpack_data(void* buffer, int bufsize, int* position) {
            Vec acorr;
            float comp_cost;
            MyMPI_Unpack(buffer, bufsize, position, &acorr[0], ndim_,
                         MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &comp_cost, 1, MPI_FLOAT);
            _p->set_gravitational_acceleration(
                    _p->get_gravitational_acceleration() + acorr);
            _p->add_comp_cost(comp_cost);
        }
    };

    /**
     * @brief Auxiliary class used to communicate particle information to other
     * MPI processes during the Ewald correction treewalk
     */
    class Import {
      private:
        /*! @brief Position for which the Ewald correction is calculated */
        Vec _pos;

        /*! @brief Old acceleration, used for the relative opening criterion */
        double _olda;

        /*! @brief Result of the local treewalk (is exported) */
        Vec _acorr;

        /*! @brief Computational cost of the local treewalk (is exported) */
        float _comp_cost;

        /*! @brief Softening length of the external Particle */
        double _hsoft;

        /*! @brief Type of the external Particle (PARTTYPE_GAS/PARTTYPE_DM) */
        ParticleType _type;

      public:
        /**
         * @brief Constructor. Initialize the import from the given
         * communication buffer
         *
         * @param buffer Communication buffer
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        Import(void* buffer, int bufsize, int* position) {
            MyMPI_Unpack(buffer, bufsize, position, &_pos[0], ndim_,
                         MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_olda, 1, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_hsoft, 1, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_type, 1, MPI_INT);
        }

        /**
         * @brief Get the position of the external Particle
         *
         * @return Position of the external Particle
         */
        Vec get_position() {
            return _pos;
        }

        /**
         * @brief Get the old acceleration of the external Particle
         *
         * @return Old acceleration of the Particle
         */
        double get_olda() {
            return _olda;
        }

        /**
         * @brief Get the softening length of the external Particle
         *
         * @return Softening length of the external Particle
         */
        double get_hsoft() {
            return _hsoft;
        }

        /**
         * @brief Get the type of the external Particle
         *
         * @return Type of the external Particle (PARTTYPE_GAS/PARTTYPE_DM)
         */
        ParticleType get_type() {
            return _type;
        }

        /**
         * @brief Set the Ewald correction force
         *
         * @param acorr Value of the Ewald correction
         */
        void set_acorr(Vec acorr) {
            _acorr = acorr;
        }

        /**
         * @brief Set the computational cost
         *
         * @param comp_cost Value of the computational cost
         */
        void set_comp_cost(float comp_cost) {
            _comp_cost = comp_cost;
        }

        /**
         * @brief Add relevant data to the given communication stream to send
         * them back to the original process
         *
         * @param buffer Communication buffer
         * @param bufsize Size of the buffer
         * @param position Current position of the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            MyMPI_Pack(&_acorr[0], ndim_, MPI_DOUBLE, buffer, bufsize,
                       position);
            MyMPI_Pack(&_comp_cost, 1, MPI_FLOAT, buffer, bufsize, position);
        }
    };

    BHEwaldGravityWalker(Particle* p, EwaldTable& ewald_table,
                         bool local = true);
    BHEwaldGravityWalker(Import& import, EwaldTable& ewald_table);

    ~BHEwaldGravityWalker();

    Vec get_acceleration();

    bool splitnode(TreeNode* node);
    void nodeaction(TreeNode* node);
    void leafaction(Leaf* leaf);
    void pseudonodeaction(PseudoNode* pseudonode);
    bool export_to_pseudonode(PseudoNode* pseudonode);

    void set_position(Vec position);
    Vec get_position();

    bool periodicsplitnode(TreeNode* node, EwaldTable& ewald_table);
    void periodicpseudonodeaction(PseudoNode* pseudonode,
                                  EwaldTable& ewald_table);
    void periodicleafaction(Leaf* leaf, EwaldTable& ewald_table);

    void after_walk();

    void after_walk(Import& import);

    Export get_export();
};

/**
 * @brief TreeWalker specialization used to calculate gravitational potentials
 * using a Barnes-Hut opening criterion
 */
class BHPotentialWalker : public PeriodicTreeWalker {
  private:
    /*! @brief Particle for which the potential is calculated */
    Particle* _p;

    /*! @brief Position of the Particle */
    Vec _position;

    /*! @brief Current value of the gravitational potential */
    double _epot;

    /*! @brief Softening length of the Particle */
    double _hsoft;

    /*! @brief Inverse softening length of the Particle */
    double _hsoftinv;

    /*! @brief Inverse softening length of the Particle cubed */
    double _hsoftinv3;

    /*! @brief Flag indicating if the treewalk is done for a local Particle */
    bool _local;

    /*! @brief Timer for the computational cost of the treewalk */
    Timer _comp_cost;

    /*! @brief GravityKernel used for close encounter softening */
    GravityKernel* _kernel;

  public:
    /**
     * @brief Auxiliary class to communicate particle information to other MPI
     * processes during the BH gravity potential treewalk
     */
    class Export {
      private:
        /*! @brief Particle for which the gravitational potential is
         *  calculated */
        Particle* _p;

        /*! @brief Position of the Particle (is exported) */
        Vec _pos;

        /*! @brief Softening length of the Particle (is exported) */
        double _hsoft;

        /*! @brief Type of the Particle (PARTTYPE_GAS/PARTTYPE_DM) */
        ParticleType _type;

      public:
        /**
         * @brief Constructor
         *
         * @param p Particle for which the gravitational potential is calculated
         */
        Export(Particle* p) {
            _p = p;
            _pos = p->get_position();
            _hsoft = p->get_hsoft();
        }

        /**
         * @brief Add relevant data to the given communication buffer to export
         * them to another MPI-process
         *
         * @param buffer Communication buffer
         * @param bufsize Buffer size
         * @param position Current position of the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            MyMPI_Pack(&_pos[0], ndim_, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_hsoft, 1, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_type, 1, MPI_INT, buffer, bufsize, position);
        }

        /**
         * @brief Read the response from another process from the given
         * communication buffer and finalize the treewalk
         *
         * @param buffer Communication buffer
         * @param bufsize Buffer size
         * @param position Current position of the buffer (is updated)
         */
        void unpack_data(void* buffer, int bufsize, int* position) {
            double epot;
            float comp_cost;
            MyMPI_Unpack(buffer, bufsize, position, &epot, 1, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &comp_cost, 1, MPI_FLOAT);
            _p->set_gravitational_potential(_p->get_gravitational_potential() +
                                            epot);
            _p->add_comp_cost(comp_cost);
        }
    };

    /**
     * @brief Auxiliary class to communicate particle information to other MPI
     * processes during the BH gravity potential treewalk
     */
    class Import {
      private:
        /*! @brief Position of the external Particle */
        Vec _pos;

        /*! @brief Local result of the treewalk (is exported) */
        double _epot;

        /*! @brief Computational cost of the local treewalk (is exported) */
        float _comp_cost;

        /*! @brief Softening length of the external Particle */
        double _hsoft;

        /*! @brief Type of the external Particle (PARTTYPE_GAS/PARTTYPE_DM) */
        ParticleType _type;

      public:
        /**
         * @brief Constructor. Initialize the import based on the given
         * communication buffer
         *
         * @param buffer Communication buffer
         * @param bufsize Buffer size
         * @param position Current position of the buffer (is updated)
         */
        Import(void* buffer, int bufsize, int* position) {
            MyMPI_Unpack(buffer, bufsize, position, &_pos[0], ndim_,
                         MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_hsoft, 1, MPI_DOUBLE);
            MyMPI_Unpack(buffer, bufsize, position, &_type, 1, MPI_INT);
        }

        /**
         * @brief Get the position of the external Particle
         *
         * @return Position of the external Particle
         */
        Vec get_position() {
            return _pos;
        }

        /**
         * @brief Get the softening length of the external Particle
         *
         * @return Softening length of the external Particle
         */
        double get_hsoft() {
            return _hsoft;
        }

        /**
         * @brief Get the type of the external Particle
         *
         * @return Type of the external Particle (PARTTYPE_GAS/PARTTYPE_DM)
         */
        ParticleType get_type() {
            return _type;
        }

        /**
         * @brief Set the gravitational potential
         *
         * @param epot Value of the gravitational potential
         */
        void set_epot(double epot) {
            _epot = epot;
        }

        /**
         * @brief Set the computational cost
         *
         * @param comp_cost Value of the computational cost
         */
        void set_comp_cost(float comp_cost) {
            _comp_cost = comp_cost;
        }

        /**
         * @brief Add relevant data to the given communication buffer to export
         * them back to the original process
         *
         * @param buffer Communication buffer
         * @param bufsize Buffer size
         * @param position Current position of the buffer (is updated)
         */
        void pack_data(void* buffer, int bufsize, int* position) {
            MyMPI_Pack(&_epot, 1, MPI_DOUBLE, buffer, bufsize, position);
            MyMPI_Pack(&_comp_cost, 1, MPI_FLOAT, buffer, bufsize, position);
        }
    };

    BHPotentialWalker(Particle* p, bool local = true);
    BHPotentialWalker(Import& import);

    ~BHPotentialWalker();

    double get_epot();

    bool splitnode(TreeNode* node);
    void nodeaction(TreeNode* node);
    void leafaction(Leaf* leaf);
    void pseudonodeaction(PseudoNode* pseudonode);
    bool export_to_pseudonode(PseudoNode* pseudonode);

    void set_position(Vec position);
    Vec get_position();

    bool periodicsplitnode(TreeNode* node, EwaldTable& ewald_table);
    void periodicpseudonodeaction(PseudoNode* pseudonode,
                                  EwaldTable& ewald_table);
    void periodicleafaction(Leaf* leaf, EwaldTable& ewald_table);

    void after_walk();

    void after_walk(Import& import);

    Export get_export();
};

/**
 * @brief Dummy TreeWalker for analytical gravitational potentials
 *
 * The ConstantPotentialGravityWalker does nothing during the treewalk and just
 * sets the acceleration to a value given by the get_gravitational_acceleration
 * member function of the template Potential at the end of the treewalk. The
 * gravitational potential is set to the value provided by the member function
 * get_gravitational_potential.
 * This class is provided to enable the use of the same integration code as for
 * any other potential, just by switching the TreeWalker that is used.
 */
template <class T>
class ConstantPotentialGravityWalker : public PeriodicTreeWalker {
  private:
    /*! @brief Particle for which the acceleration is calculated */
    Particle* _p;

    /*! @brief Position of the Particle */
    Vec _position;

    /*! @brief Current value of the acceleration */
    Vec _a;

  public:
    /**
     * @brief Dummy Export class
     *
     * Needs to be present for general treewalks, but is never used for this
     * particular TreeWalker because there is no communication involved with
     * a fixed potential.
     */
    class Export {
      public:
        /**
         * @brief Dummy constructor
         */
        Export() {}

        /**
         * @brief Dummy export method
         *
         * @param buffer Communication buffer
         * @param bufsize Buffer size
         * @param position Position of the buffer (is not updated, because
         * nothing is done by this method)
         */
        void pack_data(void* buffer, int bufsize, int* position) {}

        /**
         * @brief Dummy import method
         *
         * @param buffer Communication buffer
         * @param bufsize Buffer size
         * @param position Current position of the buffer
         */
        void unpack_data(void* buffer, int bufsize, int* position) {}
    };

    /**
     * @brief Dummy Import class
     *
     * Needs to be present for general treewalks, but is never used for this
     * particular TreeWalker because there is no communication involved with
     * a fixed potential.
     */
    class Import {
      public:
        /**
         * @brief Dummy constructor
         *
         * @param buffer Communication buffer
         * @param bufsize Buffer size
         * @param position Current position of the buffer
         */
        Import(void* buffer, int bufsize, int* position) {}

        /**
         * @brief Dummy export method
         *
         * @param buffer Communication buffer
         * @param bufsize Buffer size
         * @param position Current position of the buffer
         */
        void pack_data(void* buffer, int bufsize, int* position) {}
    };

    /**
     * @brief Constructor
     *
     * @param p Particle for which the acceleration is calculated
     * @param local Bool indicating whether the treewalk is performed for a
     * local Particle
     */
    ConstantPotentialGravityWalker(Particle* p, bool local = true) {
        _p = p;
        _position = p->get_position();
    }

    /**
     * @brief Import constructor
     *
     * @param import Dummy import imported from another MPI-process
     */
    ConstantPotentialGravityWalker(Import& import) {
        _p = NULL;
    }

    /**
     * @brief Get the acceleration calculated by the treewalk
     *
     * @return The resulting acceleration
     */
    Vec get_acceleration() {
        return _a;
    }

    /**
     * @brief Dummy splitnode function
     *
     * @param node TreeNode on which we operate
     * @return False, since we do not walk the tree
     */
    bool splitnode(TreeNode* node) {
        return false;
    }

    /**
     * @brief Dummy nodeaction
     *
     * @param node TreeNode on which we operate
     */
    void nodeaction(TreeNode* node) {}

    /**
     * @brief Dummy leafaction
     *
     * @param leaf Leaf on which we operate
     */
    void leafaction(Leaf* leaf) {}

    /**
     * @brief Dummy pseudonodeaction
     *
     * @param pseudonode PseudoNode on which we operate
     */
    void pseudonodeaction(PseudoNode* pseudonode) {}

    /**
     * @brief Dummy export check
     *
     * @param pseudonode PseudoNode on which we operate
     * @return False, since we do not walk the tree
     */
    bool export_to_pseudonode(PseudoNode* pseudonode) {
        return false;
    }

    /**
     * @brief Set the position used for the treewalk
     *
     * @param position New position
     */
    void set_position(Vec position) {
        _position = position;
    }

    /**
     * @brief Get the position of the treewalk
     *
     * @return Position used for the treewalk
     */
    Vec get_position() {
        return _position;
    }

    /**
     * @brief Dummy periodic splitnode function
     *
     * @param node TreeNode on which we operate
     * @param ewald_table EwaldTable used for periodic correction terms
     * @return False, since we do not walk the tree
     */
    bool periodicsplitnode(TreeNode* node, EwaldTable& ewald_table) {
        return false;
    }

    /**
     * @brief Dummy periodicpseudonodeaction
     *
     * @param pseudonode PseudoNode on which we operate
     * @param ewald_table EwaldTable used for periodic correction terms
     */
    void periodicpseudonodeaction(PseudoNode* pseudonode,
                                  EwaldTable& ewald_table) {}

    /**
     * @brief Dummy periodicleafaction
     * @param leaf Leaf on which we operate
     * @param ewald_table EwaldTable used for periodic correction terms
     */
    void periodicleafaction(Leaf* leaf, EwaldTable& ewald_table) {}

    /**
     * @brief Finalize the treewalk by setting the Particle acceleration
     */
    void after_walk() {
        _a = T::get_gravitational_acceleration(_position);
        _p->set_gravitational_acceleration(_a);
        _p->set_gravitational_potential(
                T::get_gravitational_potential(_position));
    }

    /**
     * @brief Dummy walk finalizing for external treewalks
     *
     * @param import Import for which the treewalk is performed
     */
    void after_walk(Import& import) {}

    /**
     * @brief Get a dummy export for MPI communication
     *
     * @return An dummy Export that should not be used for anything, but is
     * there to make the code compile
     */
    Export get_export() {
        return Export();
    }
};

#endif
