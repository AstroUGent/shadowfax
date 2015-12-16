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
 * @file Tree.hpp
 *
 * @brief Octree used for force calculations, neighbour searches...: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TREE_HPP
#define TREE_HPP

#include <ostream>
#include <fstream>
#include <sstream>
#include "Vec.hpp"
#include "TreeRoute.hpp"
#include "EwaldTable.hpp"
#include "Hilbert.hpp"

#include "Box.hpp"
#include "Cuboid.hpp"

#include "MPIGlobal.hpp"
#include "MPIMethods.hpp"
#include "utilities/Particle.hpp"
#include "utilities/GasParticle.hpp"
#include "utilities/DMParticle.hpp"

class VorGen;

#if ndim_==3
#define numnode_ 8
#else
#define numnode_ 4
#endif

class TreeWalker;
class PeriodicTreeWalker;

/**
 * @brief Auxiliary class used to communicate TreeNode information to a
 * PseudoNode on another MPI process
 */
class NodeInfo{
private:
    /*! \brief Hilbert key of the node */
    unsigned long _key;

    /*! \brief Maximum soundspeed in the node */
    double _cmax;

    /*! \brief Maximum velocity in the node */
    double _vmax;

    /*! \brief Total mass in the node */
    double _mass;

    /*! \brief Maximum softening length in the node */
    double _hmax;

    /*! \brief Center of mass of the node */
    Vec _center_of_mass;

public:
    NodeInfo() : _key(0), _cmax(0.), _vmax(0.), _mass(0.), _hmax(0.) {}

    /**
     * @brief Constructor
     *
     * @param key Hilbert key of the node
     * @param cmax Maximum soundspeed in the node
     * @param vmax Maximum velocity in the node
     * @param mass Total mass inside the node
     * @param hmax Maximum softening length in the node
     * @param center_of_mass Center of mass of the node
     */
    NodeInfo(unsigned long key, double cmax, double vmax, double mass,
             double hmax, Vec center_of_mass) :
        _key(key), _cmax(cmax), _vmax(vmax), _mass(mass), _hmax(hmax),
        _center_of_mass(center_of_mass) {}

    /**
     * @brief Get the Hilbert key of the node
     *
     * @return Hilbert key of the node
     */
    inline unsigned long get_key(){
        return _key;
    }

    /**
     * @brief Get the maximum soundspeed in the node
     *
     * @return Maximum soundspeed in the node
     */
    inline double get_cmax(){
        return _cmax;
    }

    /**
     * @brief Get the maximum velocity in the node
     *
     * @return Maximum velocity in the node
     */
    inline double get_vmax(){
        return _vmax;
    }

    /**
     * @brief Get the total mass inside the node
     *
     * @return Total mass inside the node
     */
    inline double get_mass(){
        return _mass;
    }

    /**
     * @brief Get the maximum softening length in the node
     *
     * @return Maximum softening length in the node
     */
    inline double get_hmax(){
        return _hmax;
    }

    /**
     * @brief Get the center of mass of the node
     *
     * @return Center of mass of the node
     */
    inline Vec get_center_of_mass(){
        return _center_of_mass;
    }

    /**
     * @brief Add the node information to the given MPI buffer for communication
     * to another MPI process
     *
     * @param buffer Buffer to write to
     * @param bufsize Buffer size
     * @param position Position in the buffer (is updated)
     */
    inline void pack_data(void *buffer, int bufsize, int *position){
        MyMPI_Pack(&_key, 1, MPI_UNSIGNED_LONG, buffer, bufsize, position);
        MyMPI_Pack(&_cmax, 1, MPI_DOUBLE, buffer, bufsize, position);
        MyMPI_Pack(&_vmax, 1, MPI_DOUBLE, buffer, bufsize, position);
        MyMPI_Pack(&_mass, 1, MPI_DOUBLE, buffer, bufsize, position);
        MyMPI_Pack(&_hmax, 1, MPI_DOUBLE, buffer, bufsize, position);
        MyMPI_Pack(&_center_of_mass[0], ndim_, MPI_DOUBLE, buffer, bufsize,
                position);
    }

    /**
     * @brief MPI constructor. Initialize the node from the given communication
     * buffer
     *
     * @param buffer Buffer to read from
     * @param bufsize Buffer size
     * @param position Position in the buffer (is updated)
     */
    inline NodeInfo(void *buffer, int bufsize, int *position){
        MyMPI_Unpack(buffer, bufsize, position, &_key, 1, MPI_UNSIGNED_LONG);
        MyMPI_Unpack(buffer, bufsize, position, &_cmax, 1, MPI_DOUBLE);
        MyMPI_Unpack(buffer, bufsize, position, &_vmax, 1, MPI_DOUBLE);
        MyMPI_Unpack(buffer, bufsize, position, &_mass, 1, MPI_DOUBLE);
        MyMPI_Unpack(buffer, bufsize, position, &_hmax, 1, MPI_DOUBLE);
        MyMPI_Unpack(buffer, bufsize, position, &_center_of_mass[0], ndim_,
                MPI_DOUBLE);
    }
};

/**
 * @brief Abstract node of the Tree
 *
 * The node can be either a TreeNode, a Leaf or a PseudoNode.
 * This abstract interface defines methods that are shared by all
 * implementations and offers some utility functions to distinguish between
 * different implementations.
 */
class Node{
protected:
    /*! \brief Flag used to keep track of relevant node information */
    unsigned int _flag;

public:
    /**
     * @brief Constructor
     *
     * The first bit in the internal flag is 1 if the node is a leaf and 0
     * otherwise. The second bit is 1 if the node is the last child of its
     * parent. The second bit is 1 for a pseudonode of the tree. The third bit
     * at last is 1 for nodes that have descendants that are pseudonodes.
     *
     * @param is_leaf True if the node is a leaf of the tree, false otherwise
     * @param is_pseudo True if the node is a pseudonode of the tree
     */
    Node(bool is_leaf, bool is_pseudo = false) : _flag(is_leaf){
        _flag |= (is_pseudo<<2);
    }

    virtual ~Node(){}

    /**
     * @brief Check if the node is a leaf of the tree
     *
     * @return True if the node is a Leaf
     */
    inline bool is_leaf(){
        return (_flag&1);
    }

    /**
     * @brief Check if the node is the last child of its parent
     *
     * @return True if the node is the last child of its parent
     */
    inline bool is_last(){
        return (_flag>>1)&1;
    }

    /**
     * @brief Check if the node is a pseudonode of the tree
     *
     * @return True if the node is a PseudoNode
     */
    inline bool is_pseudo(){
        return (_flag>>2)&1;
    }

    /**
     * @brief Check if the node has only children on the local MPI process
     *
     * @return True if none of the descendants of the node is a PseudoNode
     */
    inline bool is_local(){
        return !((_flag>>3)&1);
    }

    /**
     * @brief Print the node to the given stream
     *
     * @param out std::ostream to write to
     */
    virtual void print(std::ostream& out)=0;

    /**
     * @brief Finalize the node
     *
     * A finalized Tree is reordered for maximum efficiency and all its
     * treenodes have valid properties.
     *
     * @param sibling Node on the same level that is next in the efficient Tree
     * traversal ordering
     * @param last Flag indicating whether the node is the last child of its
     * parent in the efficient Tree traversal ordering
     * @return Number of descendants of this node
     */
    virtual unsigned int finalize(Node* sibling, bool last)=0;

    /**
     * @brief Get the total mass inside this node
     *
     * @return The mass of this node
     */
    virtual double get_mass()=0;

    /**
     * @brief Get the maximum softening length in this node
     *
     * @return Maximum softening length in the node
     */
    virtual double get_hmax()=0;

    /**
     * @brief Get the given coordinate of the center of mass of the node
     *
     * @param index Index of the coordinate (x = 0, y = 1, z = 2)
     * @return Requested coordinate of the center of mass of the node
     */
    virtual double get_center_of_mass(unsigned int index)=0;

    /**
     * @brief Get the maximum soundspeed in the node
     *
     * @return Maximum soundspeed in the node
     */
    virtual double get_cmax()=0;

    /**
     * @brief Get the maximum velocity in the node
     *
     * @return Maximum velocity in the node
     */
    virtual double get_vmax()=0;

    /**
     * @brief Recursively walk the node using the given TreeWalker
     *
     * @deprecated This way of traversing the tree is no longer used
     *
     * @param walker TreeWalker used to walk the tree
     */
    virtual void walk(TreeWalker& walker)=0;
};

/**
 * @brief Pair of Node pointer and associated Hilbert key
 */
class ExportNode : public Hilbert_Object {
private:
    /*! \brief Node associated with this export node */
    Node* _node;

public:
    /**
     * @brief Constructor
     *
     * @param node Node associated with this export node
     * @param key Hilber key associated with the node
     */
    ExportNode(Node* node, unsigned long key)
        : Hilbert_Object(key), _node(node) {}

    ~ExportNode(){}

    /**
     * @brief Get the node associated with this export node
     *
     * @return The node
     */
    inline Node* get_node(){
        return _node;
    }
};

/**
 * @brief Node representing a TreeNode on another MPI process
 *
 * The Node acts as a TreeNode when it is not opened, but when it has to be
 * opened, data is communicated to another process. The properties of the
 * PseudoNode also come from another process.
 */
class PseudoNode : public Node, public Hilbert_Object {
private:
    /*! \brief Node on the same level that is next in the efficient Tree
     *  traversal order */
    Node* _sibling;

    /*! \brief Box specifying the geometrical dimensions of the node */
    Box _box;

    /*! \brief Rank of the MPI process holding the TreeNode corresponding to
     *  the PseudoNode */
    unsigned int _src;

    /*! \brief Total mass inside the node */
    double _total_mass;

    /*! \brief Center of mass of the node */
    Vec _center_of_mass;

    /*! \brief Maximum soundspeed in the node */
    double _cmax;

    /*! \brief Maximum velocity in the node */
    double _vmax;

    /*! \brief Maximum softening length in the node */
    double _hmax;

public:
    PseudoNode(Box box, unsigned int src, unsigned long key);
    ~PseudoNode();

    unsigned int finalize(Node *sibling, bool last);
    double get_mass();
    double get_hmax();
    double get_center_of_mass(unsigned int index);
    Vec get_center_of_mass_node();
    double get_distance(Vec& coords);
    double get_cmax();
    double get_vmax();
    void print(std::ostream& out);
    unsigned int get_source();

    void set_mass(double mass);
    void set_hmax(double hmax);
    void set_center_of_mass(Vec com);
    void set_cmax(double cmax);
    void set_vmax(double vmax);

    double get_width();
    Vec get_center();

    Node* get_sibling();

    void walk(TreeWalker &walker);
};

/**
 * @brief Node that represents an internal node of the Tree
 *
 * The internal node has 4 or 8 children (depending on the dimension) which can
 * themselves be a Leaf or a TreeNode (or a PseudoNode). The TreeNode has
 * general properties like a center of mass and a total mass that are used
 * during treewalks when it does not have to be opened.
 */
class TreeNode : public Node{
private:
    /*! \brief Box specifying the geometrical dimensions of the node */
    Box _box;
    union{
        /*! \brief Child nodes of the treenode during tree construction */
        Node* _nodes[numnode_];

        struct{
            /*! \brief Node on the same level that is next in the efficient Tree
            *  traversal order */
            Node* _sibling;

            /*! \brief First child node of this treenode in the efficient Tree
             *  traversal order */
            Node* _child;

            /*! \brief Total mass inside this node */
            double _total_mass;

            /**
             * @brief Center of mass of the node
             *
             * This is not a Vec, since you cannot store a Vec inside a
             * nameless struct inside a union...
             */
            double _center_of_mass[ndim_];

            /*! \brief Maximum soundspeed inside the node */
            double _cmax;

            /*! \brief Maximum velocity inside the node */
            double _vmax;

            /*! \brief Mass of the node in the mesh regularization algorithm */
            double _mesh_m;

            /*! \brief Maximum softening length of the node */
            double _hmax;
        };
    };

public:
    TreeNode(Box box);
    ~TreeNode();

    void add_particle(Particle* p, unsigned int level = 0);
    void print(std::ostream& out);
    unsigned int finalize(Node* sibling, bool last);
    double get_mass();
    double get_mass_node();
    double get_hmax();
    double get_center_of_mass(unsigned int index);
    double get_center_of_mass_node(unsigned int index);
    Vec get_center_of_mass_node();
    double get_distance(Vec& coords);
    double get_cmax();
    double get_vmax();

    void set_mesh_m(double mesh_m);
    double get_mesh_m();

    Node* get_child();
    Node* get_sibling();

    double get_width();
    Vec get_center();

    bool interval_overlap(unsigned long a1, unsigned long a2, unsigned long b1,
                          unsigned long b2);
    void add_pseudoparticles(unsigned long keylow, unsigned long keyhigh,
                             unsigned int pindex, unsigned int level,
                             unsigned long key,
                             std::vector<PseudoNode*>& pseudonodes);
    void get_exportnodes(unsigned long keylow, unsigned long keyhigh,
                         unsigned int level, unsigned long key,
                         std::vector<ExportNode*>& exportnodes);
    void update_quantities(bool all = false);
    void update_mesh_masses();

    void walk(TreeWalker &walker);
};

/**
 * @brief Node that represents a Leaf of the Tree
 *
 * The Leaf stores a pointer to a Particle and provides functions to access
 * Particle properties.
 */
class Leaf : public Node{
private:
    /*! \brief Particle stored in this leaf */
    Particle* _particle;

    /*! \brief Node on the same level that is next in the efficient Tree
     *  traversal order */
    Node* _sibling;

public:
    Leaf();
    ~Leaf();

    void add_particle(Particle* p);
    Particle* get_particle();
    void print(std::ostream& out);
    unsigned int finalize(Node* sibling, bool last);
    double get_mass();
    double get_hmax();
    double get_center_of_mass(unsigned int index);
    double get_vmax();
    double get_cmax();

    Node* get_sibling();

    void walk(TreeWalker &walker);
};

/**
 * @brief Particle octree
 *
 * Used for neighbour finding, hydrodynamical timestep calculation and
 * gravitational treewalks.
 * The Tree is constructed using a top-down incremental construction algorithm
 * that makes extensive use of the Hilbert ordering of the particles. For
 * parallel versions of the program, it is split over several processes by
 * making use of pseudonodes (Springel 2005), each containing a fraction of
 * the total space filling Hilbert curve.
 */
class Tree{
private:
    /*! \brief Ewald table used to calculate periodic force corrections */
    EwaldTable *_ewald_table;

    /*! \brief Root node of the tree */
    Node* _root;

    /*! \brief Box specifying the geometric dimensions of the tree */
    Cuboid _box;

    /*! \brief Side length of the tree */
    double _side;

    /*! \brief Hilbert key range of nodes that are stored on the local MPI
     *  process */
    unsigned long _keyrange[2];

    /*! \brief List of pseudonodes that are stored on other MPI processes */
    std::vector<PseudoNode*> _pseudonodes;

    /*! \brief List of nodes that have pseudonodes on other MPI processes */
    std::vector<ExportNode*> _exportnodes;

    /*! \brief Flag indicating if the simulation box is periodic (true) or
     *  reflective (false) */
    bool _periodic;

public:
    Tree(Cuboid box, bool periodic = false, bool do_ewald = false,
         double alpha = 2., unsigned int size = 64);
    ~Tree();

    void set_periodic(bool periodic);

    void add_particle(Particle* p);
    void print(std::ostream& out);
    void print_top_quantities(std::ostream& out);
    unsigned int finalize();

    void get_neighbours(Vec& coords, double radius,
                        std::vector<GasParticle*>& ngblist,
                        std::vector<bool>& exportlist);
    void get_neighbours_outside(Vec& coords, double radius,
                                std::vector<VorGen*>& ngblist,
                                std::vector<bool>& exportlist,
                                unsigned int index=0);

    Particle* get_closest(Vec& coords, double r);
    Particle* get_periodic_closest(Vec& coords, double r);

    void set_velocities();
    void set_mesh_masses();

    /**
     * @brief Get the list of pseudonodes in the local tree
     *
     * @return PseudoNodes in the local tree
     */
    std::vector<PseudoNode> get_pseudonode_info(){
        std::vector<PseudoNode> pseudos;
        return pseudos;
    }

    /**
     * @brief Set the info corresponding to the given pseudonodes
     *
     * @deprecated Method does nothing and is not used.
     *
     * @param pseudonodes List of PseudoNodes
     */
    void set_pseudonode_info(std::vector<PseudoNode>& pseudonodes){}

    void reset(Cuboid box);

    void set_keyrange(unsigned long keylow, unsigned long keyhigh);
    void add_pseudoparticles(unsigned long keylow, unsigned long keyhigh,
                             unsigned int pindex);

    void exchange_pseudonodes();

    /**
     * @brief Walk the tree efficiently using the given TreeWalker
     *
     * @param walker TreeWalker containing the necessary functions to walk the
     * tree
     */
    template<typename T> void walk_tree(T& walker){
        if(_periodic){
            walker.set_boxsize(_side);
        }
        Node* stack = ((TreeNode*)_root)->get_child();
        while(stack){
            if(stack->is_leaf()){
                if(stack->is_pseudo()){
                    walker.pseudonodeaction(((PseudoNode*)stack));
                    stack = ((PseudoNode*)stack)->get_sibling();
                } else {
                    walker.leafaction((Leaf*)stack);
                    stack = ((Leaf*)stack)->get_sibling();
                }
            } else {
                if(walker.splitnode((TreeNode*)stack)){
                    stack = ((TreeNode*)stack)->get_child();
                } else {
                    stack = ((TreeNode*)stack)->get_sibling();
                }
            }
        }
    }

    /**
     * @brief Walk the tree efficiently using the given TreeWalker and keep
     * track of MPI communications
     *
     * @param walker TreeWalker containing the necessary functions to walk the
     * tree
     * @param exports List of Export objects that have to be communicated to
     * other MPI processes to perform external parts of the treewalk
     */
    template<typename T> void walk_tree(
            T& walker, std::vector< std::vector<typename T::Export> >& exports
            ){
        if(_periodic){
            walker.set_boxsize(_side);
        }
        // we cannot just add Exports to exports, since we only want to add
        // every Export at most once for every process
        std::vector<bool> do_export(exports.size(), false);
        Node* stack = ((TreeNode*)_root)->get_child();
        while(stack){
            if(stack->is_leaf()){
                if(stack->is_pseudo()){
                    if(walker.export_to_pseudonode((PseudoNode*)stack)){
                        do_export[((PseudoNode*)stack)->get_source()] = true;
                    }
                    stack = ((PseudoNode*)stack)->get_sibling();
                } else {
                    walker.leafaction((Leaf*)stack);
                    stack = ((Leaf*)stack)->get_sibling();
                }
            } else {
                if(walker.splitnode((TreeNode*)stack)){
                    stack = ((TreeNode*)stack)->get_child();
                } else {
                    stack = ((TreeNode*)stack)->get_sibling();
                }
            }
        }
        for(unsigned int i = 0; i < exports.size(); i++){
            if(do_export[i]){
                exports[i].push_back(walker.get_export());
            }
        }
    }

    /**
     * @brief Same as Tree::walk_tree(), but for periodic treewalks
     *
     * @deprecated This method possibly no longer works
     *
     * @param walker TreeWalker containing the necessary functions to walk the
     * tree
     */
    template<typename T> void walk_tree_periodic(T& walker){
        walk_tree(walker);
        // now do the periodic stuff
        Node* stack = ((TreeNode*)_root)->get_child();
        while(stack){
            if(stack->is_leaf()){
                if(stack->is_pseudo()){
                    walker.periodicpseudonodeaction(((PseudoNode*)stack),
                                                    _ewald_table);
                    stack = ((PseudoNode*)stack)->get_sibling();
                } else {
                    walker.periodicleafaction((Leaf*)stack, _ewald_table);
                    stack = ((Leaf*)stack)->get_sibling();
                }
            } else {
                if(walker.periodicsplitnode((TreeNode*)stack, _ewald_table)){
                    stack = ((TreeNode*)stack)->get_child();
                } else {
                    stack = ((TreeNode*)stack)->get_sibling();
                }
            }
        }
    }

    /*! \brief Tags used to distinguish between MPI messages */
    enum MPI_TAGS{
        /*! Tag for the sending of Export information to external processes */
        TAG_EXPORT,
        /*! Tag for the receiving of Import information on the local process */
        TAG_IMPORT
    };

    /**
     * @brief Generic function to walk the tree with a given TreeWalker on a
     * given set of particles
     *
     * @param particles ParticleVector containing a set of particles
     * @param gas True if the treewalk should be done for the GasParticles in
     * the ParticleVector
     * @param dm True if the treewalk should be done fot the DMParticles in the
     * ParticleVector
     * @param current_time unsigned long integer current simulation time
     */
    template<class T, typename P> void walk_tree(P& particles, bool gas,
                                                 bool dm,
                                                 unsigned long current_time){
        // the walker has two associated classes called Import and Export, which
        // are used to communicate relevant information to and from other
        // processes
        typedef typename T::Export TExport;
        typedef typename T::Import TImport;
        std::vector< std::vector<TExport> > exports(MPIGlobal::size);
        if(gas){
            for(unsigned int i = 0; i < particles.gassize(); i++){
                if(particles.gas(i)->get_endtime() == current_time){
                    T walker(particles.gas(i));
                    walk_tree(walker, exports);
                    walker.after_walk();
                }
            }
        }

        if(dm){
            // first do the local treewalk. If the treewalk encounters a part
            // of the three that is on another process, an Export is added to
            // the list
            for(unsigned int i = 0; i < particles.dmsize(); i++){
                if(particles.dm(i)->get_endtime() == current_time){
                    T walker(particles.dm(i));
                    walk_tree(walker, exports);
                    walker.after_walk();
                }
            }
        }

        // now do the parts of the treewalk on other processes
        // for this, we first communicate the TExports from the local treewalk
        // to the respective processes. For every TExport that is received, a
        // TImport is created on the local process and used to do a treewalk.
        // The TImports are then sent back to their original processes. Finally,
        // the received TImports are used to update the local TExports,
        // which takes into account the contribution of other processes to the
        // treewalk.
        // All communications are non-blocking. Furthermore, we use MPI_Probe to
        // query the status of messages. This way, we do not have to communicate
        // the sizes of the buffers explicitly: we just use the size provided by
        // the MPI_Status object.

        // Note to self: in this case, it might be a good idea to use a non
        // blocking MPI_IProbe. We now have to check manually if a send has
        // completed, but it might be easier to just launch send and receives
        // and wait for any of them to finish
        if(MPIGlobal::size > 1){
            // MPI_Requests used for the non-blocking sends. They are never
            // actually used, but we need to provide them for calls to MPI_Isend
            // And we need to free them with an MPI_Waitall or MPI_Test to make
            // sure the associated memory is not leaked!!
            std::vector<MPI_Request> reqs((MPIGlobal::size-1)*2,
                                          MPI_REQUEST_NULL);
            // Offsets in the MPIGlobal buffers. For every external process, we
            // reserve 2 blocks in each buffer to allow for non-blocking
            // communication of buffers (if we would use a single buffer, then
            // subsequent sends would overwrite it)
            std::vector<unsigned int> buffers((MPIGlobal::size-1)*2);
            unsigned int bufsize = MPIGlobal::sendsize/buffers.size();
            for(unsigned int i = 0; i < buffers.size(); i++){
                buffers[i] = i*bufsize;
            }

            unsigned int typesize = std::max(sizeof(TExport), sizeof(TImport));
            unsigned int maxsize = bufsize/typesize;
            if(!(bufsize%typesize)){
                maxsize--;
            }

            // pack the local TExports and send them to their respective
            // processes
            vector<int> allflag(MPIGlobal::size-1);
            vector<unsigned int> numsend(MPIGlobal::size-1);
            vector<unsigned int> sendsizes(MPIGlobal::size-1);
            // we keep track of the total number of messages that should be sent
            // and received. We send and receive at least one message for every
            // process. For every message that is sent, an answer has to be
            // received. For incomplete sends, we need to send additional
            // messages. For incomplete receives, we have to receive additional
            // messages.
            int numtoreceive = MPIGlobal::size-1;
            int numtosend = MPIGlobal::size-1;
            int numrecv = 0;
            int numsent = 0;
            for(int i = 0; i < MPIGlobal::size-1; i++){
                sendsizes[i] =
                        exports[(MPIGlobal::rank+1+i)%MPIGlobal::size].size();

                int send_pos = 0;
                for(unsigned int si = 0; si < std::min(maxsize, sendsizes[i]);
                    si++){
                    exports[(MPIGlobal::rank+1+i)%MPIGlobal::size][si].
                            pack_data(&MPIGlobal::sendbuffer[buffers[2*i]],
                                      bufsize, &send_pos);
                }
                allflag[i] = (sendsizes[i] <= maxsize);
                if(!allflag[i]){
                    numtosend++;
                }
                numsend[i] = 1;
                // add continuation signal
                MyMPI_Pack(&allflag[i], 1, MPI_INT,
                           &MPIGlobal::sendbuffer[buffers[2*i]], bufsize,
                           &send_pos);

                MyMPI_Isend(&MPIGlobal::sendbuffer[buffers[2*i]], send_pos,
                            MPI_PACKED, (MPIGlobal::rank+1+i)%MPIGlobal::size,
                            MPI_TAGS::TAG_EXPORT, &reqs[2*i]);
                numsent++;
                numtoreceive++;
            }

            // receive loop: receive 2 messages from every external process: a
            // buffer with TExports and a buffer with TImports
            // We use MPI_Probe to detect arriving messages, independent of
            // their origin or tag. We then receive and process these messages
            // according to their tag.
            vector<unsigned int> numreceived(MPIGlobal::size-1, 0);
            while(numrecv < numtoreceive || numsent < numtosend){
                MPI_Status status;
                if(numsent < numtosend){
                    for(int j = 0; j < MPIGlobal::size-1; j++){
                        if(!allflag[j]){
                            int flag;
                            MyMPI_Test(&reqs[j], &flag, &status);
                            if(flag){
                                int send_pos = 0;
                                for(unsigned int si = numsend[j]*maxsize;
                                    si < std::min((numsend[j]+1)*maxsize,
                                                  sendsizes[j]);
                                    si++){
                                    exports[(MPIGlobal::rank+1+j)%
                                            MPIGlobal::size][si].pack_data(
                                                &MPIGlobal::sendbuffer
                                                [buffers[2*j]], bufsize,
                                                &send_pos
                                                );
                                }
                                allflag[j] = (sendsizes[j] <=
                                              (numsend[j]+1)*maxsize);
                                if(!allflag[j]){
                                    numtosend++;
                                }
                                numsend[j]++;
                                // add continuation signal
                                MyMPI_Pack(&allflag[j], 1, MPI_INT,
                                           &MPIGlobal::sendbuffer[buffers[2*j]],
                                           bufsize, &send_pos);

                                MyMPI_Isend(&MPIGlobal::sendbuffer
                                            [buffers[2*j]], send_pos,
                                            MPI_PACKED, (MPIGlobal::rank+1+j)%
                                            MPIGlobal::size,
                                            MPI_TAGS::TAG_EXPORT, &reqs[2*j]);
                                numsent++;
                                numtoreceive++;
                            }
                        }
                    }
                }

                if(numrecv < numtoreceive){
                    int index;
                    int tag;
                    int recv_pos;
                    int send_pos;

                    // wait for a message from any source and with any tag to
                    // arrive
                    MyMPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, &status);
                    numrecv++;
                    // query the MPI_Status object to retrieve tag, source and
                    // size
                    tag = status.MPI_TAG;
                    index = status.MPI_SOURCE;
                    int nelements;
                    MyMPI_Get_count(&status, MPI_PACKED, &nelements);
                    // select an index in the buffers to use for receiving and
                    // sending
                    int freebuffer;
                    if(index == MPIGlobal::size-1){
                        freebuffer = 2*MPIGlobal::rank;
                    } else {
                        freebuffer = 2*index;
                    }
                    // TAG_EXPORT: the arriving message contains TExports. We
                    // receive them as TImports and walk the local tree. We then
                    // send the TImports back to their original process
                    if(tag == MPI_TAGS::TAG_EXPORT){
                        MyMPI_Recv(&MPIGlobal::recvbuffer[buffers[freebuffer]],
                                   nelements, MPI_PACKED, index,
                                   MPI_TAGS::TAG_EXPORT, &status);
                        recv_pos = 0;
                        send_pos = 0;
                        while(recv_pos < nelements-4){
                            TImport import(
                                        &MPIGlobal::recvbuffer
                                        [buffers[freebuffer]], nelements,
                                        &recv_pos
                                    );
                            T walker(import);
                            particles.get_tree().walk_tree(walker);
                            walker.after_walk(import);
                            import.pack_data(&MPIGlobal::sendbuffer
                                             [buffers[freebuffer+1]], bufsize,
                                             &send_pos);
                        }
                        int flag;
                        MyMPI_Unpack(&MPIGlobal::recvbuffer
                                     [buffers[freebuffer]], nelements,
                                     &recv_pos, &flag, 1, MPI_INT);
                        if(!flag){
                            numtoreceive++;
                        }

                        MyMPI_Isend(&MPIGlobal::sendbuffer
                                    [buffers[freebuffer+1]], send_pos,
                                    MPI_PACKED, index, MPI_TAGS::TAG_IMPORT,
                                    &reqs[freebuffer+1]);
                    }

                    // TAG_IMPORT: the arriving message contains TImports. We
                    // unpack them in their respective TExport to add the
                    // contribution of the source process to the local particles
                    if(tag == MPI_TAGS::TAG_IMPORT){
                        MyMPI_Recv(&MPIGlobal::recvbuffer
                                   [buffers[freebuffer+1]], nelements,
                                   MPI_PACKED, index, MPI_TAGS::TAG_IMPORT,
                                   &status);
                        unsigned int j = numreceived[freebuffer/2];
                        recv_pos = 0;
                        while(recv_pos < nelements){
                            exports[index][j].unpack_data(
                                        &MPIGlobal::recvbuffer
                                        [buffers[freebuffer+1]], nelements,
                                        &recv_pos
                                    );
                            j++;
                        }
                        numreceived[freebuffer/2] = j;
                    }
                }
            }
            // We do not have to do an extra MPI_Test loop, since we know that
            // every original send corresponds to a receive above

            // Do this to prevent communications from overlapping
            // If process 0 finishes first and enters another treewalk
            // and again sends to process 1 e.g., process 1 might receive
            // this message instead of the message of process 2 it was waiting
            // for. Result: process 1 will try to unpack data that are not
            // packed as expected and will crash
            // By putting a barrier here, we prevent this from happening
            // Another solution would be to use treewalk specific message tags,
            // but this is difficult when using template walkers.
            vector<MPI_Status> status((MPIGlobal::size-1)*2);
            MyMPI_Waitall((MPIGlobal::size-1)*2, &reqs[0], &status[0]);
            MyMPI_Barrier();
        }
    }
};

#endif // TREE_HPP
