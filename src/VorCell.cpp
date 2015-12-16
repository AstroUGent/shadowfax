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
 * @file VorCell.cpp
 *
 * @brief Voronoi cell for old algorithm: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "VorCell.hpp"
#include "Simplex.hpp"
#include "ExArith.h"
#include "VorFace.hpp"
#include "utilities/GasParticle.hpp"
#include "Error.hpp"
#include <algorithm>
#include <vector>
#include <set>
#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace std;

/**
 * @brief Initialize VoronoiCell around the given central point
 *
 * @param point The VorGen which is the center of the cell
 */
VorCell::VorCell(VorGen* point){
    _central_point = point;
    _volume = -1;
    _total_area = -1;
}

/**
 * @brief Destructor
 *
 * Does nothing.
 */
VorCell::~VorCell(){
}

/**
 * @brief Print the cell (in ascii) to the given stream
 *
 * Points are outputted as:
 * p:\\t(c1,c2,c3)\n (with the ci doubles)
 * Lines as:
 * l:\\t(c1,c2,c3)\\t(d1,d2,d3)\n
 * This method calls the print methods of the points and faces.
 *
 * @param stream std::stream to write to
 */
void VorCell::print(ostream& stream){
    for(unsigned int i = 0; i < _faces.size(); i++){
        _central_point->print(stream);
        stream << "\n";
        _faces[i]->print(stream);
    }
}

/**
 * @brief Print the cell to the given stream in a format that can be plotted
 * with gnuplot
 *
 * @param stream std::ostream to write to
 * @param id Index of the cell in the VorTess cell list
 */
void VorCell::print_gnuplot(ostream& stream, unsigned int id){
    _central_point->print(stream);
    stream << "\t" << id << "\n\n";
    for(unsigned int i = 0; i < _faces.size(); i++){
        _faces[i]->print(stream);
    }
    stream << "\n\n";
}

/**
 * @brief Print the cell to the given stream in a format that can be plotted
 * using POVRAY
 *
 * @param stream std::ostream to write to
 */
void VorCell::print_pov(ostream &stream){
#if ndim_==3
    stream << "sphere{<" << _central_point->x() << "," << _central_point->y()
           << "," << _central_point->z() << ">,R\n texture { Generator } }"
           << endl;
    for(unsigned int i = 0; i < _faces.size(); i++){
        _faces[i]->print_pov_frame(stream);
    }
    stream << "mesh {" << endl;
    for(unsigned int i = 0; i < _faces.size(); i++){
        _faces[i]->print_pov(stream);
    }
    stream << "}" << endl;
#endif
}

#ifndef ICMAKER
/**
 * @brief Print the cell to the given stream in a format that can be plotted
 * with Leaflet
 *
 * @param vstream std::ostream to write to
 * @param ox x-coordinate of the origin of the screen coordinate system
 * @param oy y-coordinate of the origin of the screen coordinate system
 * @param colormap ColorMap used to convert values to colors
 * @param maxW Maximal primitive variables
 * @param minW Minimal primitive variables
 */
void VorCell::print_leaflet(ostream &vstream, int ox, int oy,
                            ColorMap* colormap, StateVector maxW,
                            StateVector minW){
    short numfaces = (short)_faces.size();
    vstream.write((char*)&numfaces, sizeof(short));
    for(unsigned int i = 0; i < _faces.size(); i++){
        _faces[i]->print_leaflet(vstream, ox, oy, _central_point);
    }
    // write density
    int icolors[3];
    colormap->get_color(
                (_central_point->get_particle()->get_Wvec().rho()-
                 minW.rho())/(maxW.rho()-minW.rho()), icolors
                );
    short colors[3];
    for(int i = 0; i < 3; i++){
        colors[i] = (short)icolors[i];
    }
    vstream.write((char*)&colors[0], 3*sizeof(short));
    // write pressure
    colormap->get_color(
                (_central_point->get_particle()->get_Wvec().p()-
                 minW.p())/(maxW.p()-minW.p()), icolors
                );
    for(int i = 0; i < 3; i++){
        colors[i] = (short)icolors[i];
    }
    vstream.write((char*)&colors[0], 3*sizeof(short));
}
#endif

/**
 * @brief Check if one of the faces of the cell overlaps with the given box
 *
 * @param box Coordinates of a box
 * @return True if a face of the cell overlaps with the box
 */
bool VorCell::overlap(double *box){
    for(unsigned int i = 0; i < _faces.size(); i++){
        if(_faces[i]->overlap(box)){
           return true;
        }
    }
    return false;
}

/**
 * @brief Print the cell to the given stream in VTK-format
 *
 * @param vstream std::ostream to write vertices to
 * @param numv Number of vertices written to the vertices stream (updated)
 * @param pstream std::ostream to write polygons to
 * @param nump Number of polygons written to the polygon stream (updated)
 * @param numc Number of connections written to the polygon stream (updated)
 * @param dstream std::ostream to write densities to
 */
void VorCell::print_vtk(ostream &vstream, unsigned int &numv, ostream &pstream,
                        unsigned int &nump, unsigned int &numc,
                        ostream &dstream){
#if ndim_==2
    map<double, VorGen*> vertices;
    for(unsigned int i = _faces.size(); i--;){
        vector<VorGen*> face = _faces[i]->get_vertices();
        for(unsigned int j = 2; j--;){
            double angle = atan2(face[j]->y()-_central_point->y(),
                                 face[j]->x()-_central_point->x());
            vertices[angle] = face[j];
        }
    }
    unsigned int vsize = vertices.size();
    pstream << vsize;
    unsigned int i = 0;
    for(map<double, VorGen*>::iterator it = vertices.begin();
        it != vertices.end(); it++){
        vstream << it->second->x() << " " << it->second->y() << " 0\n";
        pstream << " " << numv+(i++);
    }
    numv += vertices.size();
    nump++;
    numc += vsize + 1;
    pstream << "\n";
    dstream << _central_point->get_particle()->get_Wvec().rho() << "\n";
#endif
}

/**
 * @brief Get the centroid of the cell
 *
 * The centroid should be set or calculated first, using VorCell::set_centroid()
 * or VorCell::calculate_centroid().
 *
 * @return Reference to the centroid of the cell
 */
Vec& VorCell::get_centroid(){
    return _centroid;
}

/**
 * @brief Set the centroid of the cell
 *
 * This version sets the centroid using an array.
 *
 * @param centroid Components of the new centroid of the cell
 */
void VorCell::set_centroid(double* centroid){
#if ndim_==3
    _centroid.set(centroid[0], centroid[1], centroid[2]);
#else
    _centroid.set(centroid[0], centroid[1]);
#endif
}

/**
 * @brief Set the centroid of the cell
 *
 * @param centroid New centroid of the cell
 */
void VorCell::set_centroid(Vec centroid){
    _centroid = centroid;
}

/**
 * @brief Calculate the centroid of the cell
 *
 * The centroid is given by a weighted mean value of the centroids of tetrahedra
 * formed by the central point and three neighbouring vertices. These are
 * weighted by their respective volumes.
 */
void VorCell::calculate_centroid(){
#if ndim_==3
    double C_tot[3], A_tot;
    C_tot[0] = 0;
    C_tot[1] = 0;
    C_tot[2] = 0;
    A_tot = 0;
    for(unsigned int i = 0; i < _faces.size(); i++){
        vector<VorGen*> points = _faces[i]->get_vertices();
        for(unsigned int j = 1; j < points.size()-1; j++){
            double C[3], A;
            C[0] = 0.25*(points[0]->x() + points[j]->x() + points[j+1]->x() +
                    _central_point->x());
            C[1] = 0.25*(points[0]->y() + points[j]->y() + points[j+1]->y() +
                    _central_point->y());
            C[2] = 0.25*(points[0]->z() + points[j]->z() + points[j+1]->z() +
                    _central_point->z());
            double r1[3], r2[3], r3[3];
            r1[0] = points[0]->x() - points[j]->x();
            r1[1] = points[0]->y() - points[j]->y();
            r1[2] = points[0]->z() - points[j]->z();
            r2[0] = points[j]->x() - points[j+1]->x();
            r2[1] = points[j]->y() - points[j+1]->y();
            r2[2] = points[j]->z() - points[j+1]->z();
            r3[0] = points[j+1]->x() - _central_point->x();
            r3[1] = points[j+1]->y() - _central_point->y();
            r3[2] = points[j+1]->z() - _central_point->z();
            A = fabs(r1[0]*r2[1]*r3[2] + r1[1]*r2[2]*r3[0] + r1[2]*r2[0]*r3[1] -
                    r1[2]*r2[1]*r3[0] - r2[2]*r3[1]*r1[0] - r3[2]*r1[1]*r2[0]);
            A /= 6;
            C_tot[0] += C[0]*A;
            C_tot[1] += C[1]*A;
            C_tot[2] += C[2]*A;
            A_tot += A;
        }
    }
    C_tot[0] /= A_tot;
    C_tot[1] /= A_tot;
    C_tot[2] /= A_tot;
    _centroid.set(C_tot[0], C_tot[1], C_tot[2]);
#else
    double C_tot[2], A_tot;
    C_tot[0] = 0;
    C_tot[1] = 0;
    A_tot = 0;
    for(unsigned int i = 0; i < _faces.size(); i++){
        vector<VorGen*> points = _faces[i]->get_vertices();
        double C[2];
        C[0] = (points[0]->x() + points[1]->x() + _central_point->x())/3.;
        C[1] = (points[0]->y() + points[1]->y() + _central_point->y())/3.;
        double a[2], b[2];
        a[0] = points[0]->x() - _central_point->x();
        a[1] = points[0]->y() - _central_point->y();
        b[0] = points[1]->x() - _central_point->x();
        b[1] = points[1]->y() - _central_point->y();
        double A = 0.5*fabs(a[0]*b[1] - a[1]*b[0]);
        C_tot[0] += A*C[0];
        C_tot[1] += A*C[1];
        A_tot += A;
    }
    C_tot[0] /= A_tot;
    C_tot[1] /= A_tot;
    _centroid.set(C_tot[0], C_tot[1]);
#endif
}

/**
 * @brief Get the volume of the cell
 *
 * The volume should be set or calculated first, using VorCell::set_volume() or
 * VorCell::calculate_volume().
 *
 * @return Volume of the cell
 */
double VorCell::get_volume(){
    return _volume;
}

/**
 * @brief Set the volume of the cell
 *
 * @param volume New volume of the cell
 */
void VorCell::set_volume(double volume){
    _volume = volume;
}

/**
 * @brief Calculate the volume of the cell
 *
 * The volume is the sum of the volumes of triangles/tetrahedra making up the
 * cell.
 */
void VorCell::calculate_volume(){
#if ndim_==3
    _volume = 0;
    for(unsigned int i = 0; i < _faces.size(); i++){
        vector<VorGen*> points = _faces[i]->get_vertices();
        for(unsigned int j = 1; j < points.size()-1; j++){
            double A;
            double r1[3], r2[3], r3[3];
            r1[0] = points[0]->x() - points[j]->x();
            r1[1] = points[0]->y() - points[j]->y();
            r1[2] = points[0]->z() - points[j]->z();
            r2[0] = points[j]->x() - points[j+1]->x();
            r2[1] = points[j]->y() - points[j+1]->y();
            r2[2] = points[j]->z() - points[j+1]->z();
            r3[0] = points[j+1]->x() - _central_point->x();
            r3[1] = points[j+1]->y() - _central_point->y();
            r3[2] = points[j+1]->z() - _central_point->z();
            A = fabs(r1[0]*r2[1]*r3[2] + r1[1]*r2[2]*r3[0] + r1[2]*r2[0]*r3[1] -
                    r1[2]*r2[1]*r3[0] - r2[2]*r3[1]*r1[0] - r3[2]*r1[1]*r2[0]);
            A /= 6;
            _volume += A;
        }
    }
#else
    _volume = 0.;
    for(unsigned int i = _faces.size(); i--;){
        vector<VorGen*> points = _faces[i]->get_vertices();
        double a[2], b[2];
        a[0] = points[0]->x() - _central_point->x();
        a[1] = points[0]->y() - _central_point->y();
        b[0] = points[1]->x() - _central_point->x();
        b[1] = points[1]->y() - _central_point->y();
        _volume += 0.5*fabs(a[0]*b[1] - a[1]*b[0]);
    }
#endif
}

/**
 * @brief Get the total area of all faces of the cell
 *
 * The total area is calculated once and then stored for further usage.
 *
 * @return Total area of all faces of the cell
 */
double VorCell::get_total_area(){
    if(_total_area < 0){
        calculate_total_area();
    }
    return _total_area;
}

/**
 * @brief Calculate the total area of all faces of the cell
 */
void VorCell::calculate_total_area(){
    _total_area = 0.;
    for(unsigned int i = 0; i < _faces.size(); i++){
        _total_area += _faces[i]->get_area();
    }
}

/**
 * @brief Add a neighbour point of the cell to the list
 *
 * @param ngb A VorGen that is a neighbour to this cell
 */
void VorCell::add_ngb(VorGen* ngb){
    _ngbs.push_back(ngb);
}

/**
 * @brief Access the neighbours of this cell
 *
 * @return A vector containing the neighbouring points of this cell
 */
vector<VorGen*> VorCell::get_ngbs(){
    return _ngbs;
}

/**
 * @brief Add the given face to the back of the face list
 *
 * @param face VorFace to add
 */
void VorCell::add_face(VorFace* face){
    _faces.push_back(face);
}

/**
 * @brief Return the number of faces stored in this cell
 *
 * @return The number of faces stored in this cell
 */
int VorCell::number_of_faces(){
    return _faces.size();
}

/**
 * @brief Get the faces of this cell
 *
 * @return A vector containing the faces of this cell
 */
vector<VorFace*> VorCell::get_faces(){
    return _faces;
}

/**
 * @brief Estimate slope-limited gradients for the primitive variables in this
 * cell
 *
 * @param delta Arra to store the resulting gradients in
 */
void VorCell::estimate_gradient(StateVector* delta){
    StateVector W = _central_point->get_particle()->get_Wvec();
#ifndef NOMUSCL
    // gradient estimation
    StateVector Wmaxvec, Wminvec;
    Wmaxvec = W;
    Wminvec = W;
    for(unsigned int j = _faces.size(); j--;){
      double Aij = _faces[j]->get_area();
      if(Aij){
        Vec& midface = _faces[j]->get_midpoint();
        Vec c = midface - 0.5*(_central_point->get_position() +
                               _ngbs[j]->get_position());
        Vec rij = _central_point->get_position() - _ngbs[j]->get_position();
        double rnorm = rij.norm();
        StateVector Wj;
        // check if we have to mirror this particle
        if(_ngbs[j]->get_particle() == NULL){
            Wj = _central_point->get_particle()->get_Wvec();
            _faces[j]->transform(Wj, _faces[j]->get_left()->get_position(),
                                 _faces[j]->get_ngb()->get_position());
            Wj[1] = -Wj[1];
            _faces[j]->invtransform(Wj, _faces[j]->get_left()->get_position(),
                                    _faces[j]->get_ngb()->get_position());
        } else {
            Wj = _ngbs[j]->get_particle()->get_Wvec();
        }
        for(unsigned int l = ndim_; l--;){
            delta[l] += Aij*((Wj-W)*c[l]/rnorm - 0.5*(W+Wj)*rij[l]/rnorm);
        }
        Wmaxvec = max(Wmaxvec, Wj);
        Wminvec = min(Wminvec, Wj);
      }
    }
    for(unsigned int l = ndim_; l--;){
        delta[l] /= get_volume();
    }
    // slope limiting
    Vec& centroid = get_centroid();
    StateVector alphavec(1.);
    for(unsigned int l = _faces.size(); l--;){
      if(_faces[l]->get_area()){
        Vec& midface = _faces[l]->get_midpoint();
        Vec d = midface - centroid;
        StateVector deltap = delta[0]*d[0] + delta[1]*d[1];
#if ndim_==3
        deltap += delta[2]*d[2];
#endif
        alphavec = min(alphavec, max((Wmaxvec-W)/deltap, (Wminvec-W)/deltap));
      }
    }
    for(unsigned int m = ndim_; m--;){
        delta[m] *= alphavec;
    }

#else
    for(unsigned int i = 0; i < ndim_; i++){
        for(unsigned int j = 0; j < ndim_+2; j++){
            delta[i][j] = 0.;
        }
    }
#endif
}

/**
 * @brief Get the face corresponding to the given neighbour
 *
 * If the given point is not a neighbour of this cell, the code will abort.
 *
 * @param point VorGen that is a valid neighbour of this cell
 * @return VorFace corresponding to the given neighbour
 */
VorFace* VorCell::get_face(VorGen* point){
    unsigned int i = _ngbs.size();
    while(i-- && _ngbs[i] != point){}
    if(i && _ngbs[i] != point){
        cout << "Trying to ask face for a ngb that is not a neighbour!" << endl;
        my_exit();
    }
    return _faces[i];
}

/**
 * @brief Add the given index to the back of the neighbour index list
 *
 * @param ngb Index of a VorGen in the DelTess point list
 */
void VorCell::add_ngb_id(unsigned int ngb){
    _ngb_ids.push_back(ngb);
}

/**
 * @brief Get the indices of the neighbours of this cell
 *
 * @return Vector with indices of points in the DelTess point list
 */
vector<unsigned int> VorCell::get_ngb_ids(){
    return _ngb_ids;
}

/**
 * @brief Add the given index to the back of the face index list
 *
 * @param face Index of a VorFace in the VorTess face list
 */
void VorCell::add_face_id(unsigned int face){
    _face_ids.push_back(face);
}

/**
 * @brief Get the indices of the faces of this cell
 *
 * @return Vector with indices of faces in the VorTess face list
 */
vector<unsigned int> VorCell::get_face_ids(){
    return _face_ids;
}

#if ndim_==3
/**
 * @brief Calculate the volume of overlap between this cell an the given cube
 *
 * @warning Not implemented!
 *
 * @param corner Corner of the cube
 * @param side Side length of the cube
 * @return Volume of the intersection of this cell with the given cube
 */
double VorCell::overlap_volume(double *corner, double side){
    // not implemented yet
    return 0.;
}
#else
/**
 * @brief Calculate the area of overlap between this cell an the given square
 *
 * @param corner Corner of the square
 * @param side Side length of the square
 * @return Volume of the intersection of this cell with the given square
 */
double VorCell::overlap_volume(double *corner, double side){
    vector< vector<double> > faces;
    for(unsigned int i = _faces.size(); i--;){
        vector<double> points;
        vector<VorGen*> vertices = _faces[i]->get_vertices();
        points.push_back(vertices[0]->x());
        points.push_back(vertices[0]->y());
        points.push_back(vertices[1]->x());
        points.push_back(vertices[1]->y());
        faces.push_back(points);
    }

    // the bottom line of the square
    vector<double> xbottom;
    for(unsigned int i = faces.size(); i--;){
        bool a = faces[i][1] >= corner[1];
        bool b = faces[i][3] >= corner[1];
        if(a){
            if(!b){
                // the case faces[i][1] == faces[i][3] is impossible, because
                // then b should be true
                double x = (faces[i][2]-faces[i][0])/(faces[i][3]-faces[i][1])*
                        (corner[1]-faces[i][1]) + faces[i][0];
                faces[i][2] = x;
                faces[i][3] = corner[1];
                xbottom.push_back(x);
                xbottom.push_back(corner[1]);
            }
        } else {
            if(b){
                // the case faces[i][1] == faces[i][3] is impossible in this
                // case too (because then a should be true)
                double x = (faces[i][2]-faces[i][0])/(faces[i][3]-faces[i][1])*
                        (corner[1]-faces[i][1]) + faces[i][0];
                faces[i][0] = x;
                faces[i][1] = corner[1];
                xbottom.push_back(x);
                xbottom.push_back(corner[1]);
            } else {
                // the face is completely underneath the line: erase it
                faces.erase(faces.begin()+i);
            }
        }
    }
    // if one of the faces happened to be on the line, xbottom is empty
    if(xbottom.size()){
        if(xbottom.size() > 4){
            cout << "Something strange happened while determining the overlap "
                    "between a grid cell and a voronoi cell:" << endl;
            cout << "There are more than 2 intersections between the cell and "
                    "a line..." << endl;
            my_exit();
        }
        faces.push_back(xbottom);
    }

    // the top line of the square
    vector<double> xtop;
    for(unsigned int i = faces.size(); i--;){
        bool a = faces[i][1] <= corner[1]+side;
        bool b = faces[i][3] <= corner[1]+side;
        if(a){
            if(!b){
                // the case faces[i][1] == faces[i][3] is impossible, because
                // then b should be true
                double x = (faces[i][2]-faces[i][0])/(faces[i][3]-faces[i][1])*
                        (corner[1]+side-faces[i][1]) + faces[i][0];
                faces[i][2] = x;
                faces[i][3] = corner[1]+side;
                xtop.push_back(x);
                xtop.push_back(corner[1]+side);
            }
        } else {
            if(b){
                // the case faces[i][1] == faces[i][3] is impossible in this
                // case too (because then a should be true)
                double x = (faces[i][2]-faces[i][0])/(faces[i][3]-faces[i][1])*
                        (corner[1]+side-faces[i][1]) + faces[i][0];
                faces[i][0] = x;
                faces[i][1] = corner[1]+side;
                xtop.push_back(x);
                xtop.push_back(corner[1]+side);
            } else {
                // the face is completely underneath the line: erase it
                faces.erase(faces.begin()+i);
            }
        }
    }
    // if one of the faces happened to be on the line, xtop is empty
    if(xtop.size()){
        if(xtop.size() > 4){
            cout << "Something strange happened while determining the overlap "
                    "between a grid cell and a voronoi cell:" << endl;
            cout << "There are more than 2 intersections between the cell and "
                    "a line..." << endl;
            my_exit();
        }
        faces.push_back(xtop);
    }

    // the left line of the square
    vector<double> yleft;
    for(unsigned int i = faces.size(); i--;){
        bool a = faces[i][0] >= corner[0];
        bool b = faces[i][2] >= corner[0];
        if(a){
            if(!b){
                // the case faces[i][0] == faces[i][2] is impossible, because
                // then b should be true
                double y = (faces[i][3]-faces[i][1])/(faces[i][2]-faces[i][0])*
                        (corner[0]-faces[i][0]) + faces[i][1];
                faces[i][2] = corner[0];
                faces[i][3] = y;
                yleft.push_back(corner[0]);
                yleft.push_back(y);
            }
        } else {
            if(b){
                // the case faces[i][0] == faces[i][2] is impossible in this
                // case too (because then a should be true)
                double y = (faces[i][3]-faces[i][1])/(faces[i][2]-faces[i][0])*
                        (corner[0]-faces[i][0]) + faces[i][1];
                faces[i][0] = corner[0];
                faces[i][1] = y;
                yleft.push_back(corner[0]);
                yleft.push_back(y);
            } else {
                // the face is completely underneath the line: erase it
                faces.erase(faces.begin()+i);
            }
        }
    }
    // if one of the faces happened to be on the line, xtop is empty
    if(yleft.size()){
        if(yleft.size() > 4){
            cout << "Something strange happened while determining the overlap "
                    "between a grid cell and a voronoi cell:" << endl;
            cout << "There are more than 2 intersections between the cell and "
                    "a line..." << endl;
            my_exit();
        }
        faces.push_back(yleft);
    }

    // the right line of the square
    vector<double> yright;
    for(unsigned int i = faces.size(); i--;){
        bool a = faces[i][0] <= corner[0]+side;
        bool b = faces[i][2] <= corner[0]+side;
        if(a){
            if(!b){
                // the case faces[i][0] == faces[i][2] is impossible, because
                // then b should be true
                double y = (faces[i][3]-faces[i][1])/(faces[i][2]-faces[i][0])*
                        (corner[0]+side-faces[i][0]) + faces[i][1];
                faces[i][2] = corner[0]+side;
                faces[i][3] = y;
                yright.push_back(corner[0]+side);
                yright.push_back(y);
            }
        } else {
            if(b){
                // the case faces[i][0] == faces[i][2] is impossible in this
                // case too (because then a should be true)
                double y = (faces[i][3]-faces[i][1])/(faces[i][2]-faces[i][0])*
                        (corner[0]+side-faces[i][0]) + faces[i][1];
                faces[i][0] = corner[0]+side;
                faces[i][1] = y;
                yright.push_back(corner[0]+side);
                yright.push_back(y);
            } else {
                // the face is completely underneath the line: erase it
                faces.erase(faces.begin()+i);
            }
        }
    }
    // if one of the faces happened to be on the line, xtop is empty
    if(yright.size()){
        if(yright.size() > 4){
            cout << "Something strange happened while determining the overlap "
                    "between a grid cell and a voronoi cell:" << endl;
            cout << "There are more than 2 intersections between the cell and "
                    "a line..." << endl;
            my_exit();
        }
        faces.push_back(yright);
    }

    // calculate the area confined by the remaining faces
    // as the remainder of the cell should still be convex, we can use one of
    // the vertices as origin
    // the total area is then the sum of the areas of triangles formed by this
    // origin and the other faces
    double area = 0.;
    if(faces.size()){
        double origin[2] = {faces[faces.size()-1][0], faces[faces.size()-1][1]};
        for(unsigned int i = faces.size()-1; i--;){
            if(!(faces[i][0] == origin[0] && faces[i][1] == origin[1]) &&
                    !(faces[i][2] == origin[0] && faces[i][3] == origin[1])){
                // formula from wikipedia
                // (http://en.wikipedia.org/wiki/Triangle#Using_coordinates)
                area += 0.5*fabs((faces[i][0]-origin[0])*
                        (faces[i][3]-faces[i][1]) - (faces[i][0]-faces[i][2])*
                        (origin[1]-faces[i][1]));
            }
        }
    }
    return area;
}
#endif

#if ndim_==3
/**
 * @brief Same as VorCell::overlap_volume(), but for periodic boxes
 *
 * @warning Not implemented!
 *
 * @param corner Origin of the cube
 * @param side Side length of the cube
 * @return Volume of the intersection of the cell and the given cube
 */
double VorCell::periodic_overlap_volume(double *corner, double side){
    // not implemented yet
    return 0.;
}
#else
/**
 * @brief Same as VorCell::overlap_volume(), but for periodic boxes
 *
 * @warning Not implemented!
 *
 * @param corner Origin of the square
 * @param side Side length of the square
 * @return Volume of the intersection of the cell and the given square
 */
double VorCell::periodic_overlap_volume(double *corner, double side){
    double area = 0.;
    for(unsigned int i = 3; i--;){
        for(unsigned int j = 3; j--;){
            double corner_copy[2] = {corner[0]-1.+i*1., corner[1]-1.+j*1.};
            area += overlap_volume(corner_copy, side);
        }
    }
    return area;
}
#endif

/**
 * @brief Multiply eta with some cell specific local variables
 *
 * We have to do this before eta is communicated to other MPI processes, since
 * otherwise we should also communicate the cell volume to the other processes.
 */
void VorCell::finalize_eta(){
    GasParticle* p = _central_point->get_particle();
    double eta = p->get_eta()*p->get_mass()*2.8*p->get_hsoft()/6./_volume;
    p->set_eta(eta);
}

/**
 * @brief Get the gravitational force correction due to the use of variable
 * softening lengths
 *
 * @return Correction to the acceleration
 */
Vec VorCell::get_gravitational_correction(){
    Vec acorr;
    GasParticle* p = _central_point->get_particle();
    double eta_i = p->get_eta();
    for(unsigned int j = 0; j < _ngbs.size(); j++){
        GasParticle* pj;
        if( (pj = _ngbs[j]->get_particle()) ){
            double eta_j = pj->get_eta();
            Vec rij = _central_point->get_position() - _ngbs[j]->get_position();
            Vec cij = _faces[j]->get_midpoint() - 0.5*
                    (_central_point->get_position() + _ngbs[j]->get_position());
            double rijnrm = rij.norm();
            acorr += _faces[j]->get_area() * ( (eta_j - eta_i)*cij/rijnrm -
                                               0.5*(eta_i+eta_j)*rij/rijnrm );
        }
    }
    return acorr;
}

/**
 * @brief Set the characteristic length of the GasParticle corresponding to this
 * cell
 *
 * The characteristic length is the radius of a sphere/circle with the same
 * volume/face area as the cell.
 */
void VorCell::set_h(){
    double V = get_volume();
#if ndim_==3
    double h = cbrt(V*3./4./M_PI);
#else
    double h = sqrt(V/M_PI);
#endif
    _central_point->get_particle()->set_h(h);
}

/**
 * @brief Get the characteristic length of this cell
 *
 * The characteristic length is the radius of a sphere/circle with the same
 * volume/face area as the cell.
 *
 * @return Characteristic length of this cell
 */
double VorCell::get_h(){
    double V = get_volume();
#if ndim_==3
    double h = cbrt(V*3./4./M_PI);
#else
    double h = sqrt(V/M_PI);
#endif
    return h;
}

/**
 * @brief Get the velocity of the cell generator
 *
 * The velocity of the generator is equal to the local fluid velocity plus a
 * correction term that keeps the mesh regular.
 *
 * @return The velocity of the cell generator
 */
Vec VorCell::get_velocity(){
    StateVector W = _central_point->get_particle()->get_Wvec();
    Vec vel;

    if(_central_point->get_particle()->get_Wvec().rho()){
        Vec& centroid = get_centroid();
        Vec d = centroid - _central_point->get_position();
        double Ri = get_volume();
        double csnd = ((GasParticle*)_central_point->get_particle())
                ->get_soundspeed();
#if ndim_==3
        Ri = cbrt(3.*Ri/4./M_PI);
#else
        Ri = sqrt(Ri/M_PI);
#endif
        double check = 4.*d.norm()/Ri;
        if(check > 0.9){
            if(check < 1.1){
                vel = csnd/d.norm()*(d.norm()-0.9*0.25*Ri)/(0.2*0.25*Ri)*d;
            } else {
                vel = csnd/d.norm()*d;
            }
        }
    }

#if ndim_==3
    return Vec(W[1]+vel[0], W[2]+vel[1], W[3]+vel[2]);
#else
    return Vec(W[1]+vel[0], W[2]+vel[1]);
#endif
}

/**
 * @brief Get the GasParticle corresponding to this cell
 *
 * @return GasParticle corresponding to this cell
 */
GasParticle* VorCell::get_particle(){
    return _central_point->get_particle();
}

/**
 * @brief Get the triangles that make up this cell
 *
 * Used for VTK-file generation.
 *
 * @param positions std::vector with vertex positions
 * @param connectivity std::vector with connectivity information
 * @param data std::vector with StateVector data for all triangles
 */
void VorCell::get_triangles(std::vector<float> &positions,
                            std::vector<int> &connectivity,
                            std::vector<StateVector> &data){
#if ndim_==3
    for(unsigned int i = 0; i < _faces.size(); i++){
        _faces[i]->get_triangles(positions, connectivity);
        data.push_back(_central_point->get_particle()->get_Wvec());
    }
#else
    unsigned int pointsize = positions.size()/3;
    unsigned int facesize = 0;
    // vtk uses single precision instead of double precision
    // as a result, vertices for faces that are too small will overlap
    // this causes problems when plotting the vtk file in VisIt.
    // We therefore only consider faces with areas larger than 1.e-12, this
    // seems to work.
    for(unsigned int i = 0; i < _faces.size(); i++){
        if(_faces[i]->get_area() > 1.e-12){
            facesize++;
        }
    }
    connectivity.push_back(facesize);
    unsigned int lastconn = 0;
    for(unsigned int i = 0; i < _faces.size(); i++){
        if(_faces[i]->get_area() <= 1.e-12){
            continue;
        }
        vector<VorGen*> vertices = _faces[i]->get_vertices();
        // check the orientation of the face
        // if the left point for the face is the central point of this cell,
        // then the face was constructed from the point of view of this
        // cell. In this case, we use the first vertex as cell vertex. If not,
        // we use the other vertex.
        if(_faces[i]->get_left() == _central_point){
            positions.push_back(vertices[0]->x());
            positions.push_back(vertices[0]->y());
        } else {
            positions.push_back(vertices[1]->x());
            positions.push_back(vertices[1]->y());
        }
        positions.push_back(0.);
        connectivity.push_back(pointsize+lastconn);
        lastconn++;
    }
    data.push_back(_central_point->get_particle()->get_Wvec());
#endif
}
