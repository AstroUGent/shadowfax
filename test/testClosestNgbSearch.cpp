/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testClosestNgbSearch.cpp
 *
 * @brief Unit test for closest neighbour search
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "ClosestNgbSearch.hpp"
#include "myAssert.hpp"
#include "utilities/StarParticle.hpp"
#include "utilities/Tree.hpp"
using namespace std;

/**
 * @brief Unit test for closest neighbour search
 *
 * @param argc Number of command line arguments (ignored)
 * @param argv Command line arguments (ignored)
 * @return 0 on succes. Aborts otherwise
 */
int main(int argc, char** argv) {

#if ndim_ == 3
    // No intersection
    {
        Vec center(0.5, 0.5, 0.5);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.1);

        Vec anchor(0.125, 0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == false,
                  "False hit on no intersection case!");
    }
    {
        Vec center(-0.5, -0.5, -0.5);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.1);

        Vec anchor(0.125, 0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == false,
                  "False hit on no intersection case!");
    }
    // Clear intersection
    {
        Vec center(0.25, 0.3, 0.25);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.1);

        Vec anchor(0.125, 0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == true,
                  "No hit on intersection case!");
    }
    {
        Vec center(0.25, -0.05, 0.25);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.1);

        Vec anchor(0.125, 0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == true,
                  "No hit on intersection case!");
    }
    // Sphere in box
    {
        Vec center(0.2, 0.2, 0.2);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.01);

        Vec anchor(0.125, 0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == true,
                  "No hit on sphere in box case!");
    }
    // Sphere touches box corner
    {
        Vec center(0.35, 0.35, 0.35);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, sqrt(0.03));

        Vec anchor(0.125, 0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == true,
                  "No hit on sphere touches box case!");
    }
    // Sphere intersects box side, but no corner
    {
        Vec center(0.125, 0.3, 0.125);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.1);

        Vec anchor(0.125, 0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == true,
                  "No hit on sphere intersects side case!");
    }
    {
        Vec center(-0.05, 0.125, 0.125);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.1);

        Vec anchor(0.125, 0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == true,
                  "No hit on sphere intersects side case!");
    }
    // Sphere intersects periodic copy
    {
        Vec center(0.95, 0.95, 0.95);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.1);
        ngbsearch.set_boxsize(1.);

        Vec anchor(0.125, 0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == true,
                  "No hit on sphere intersects periodic copy case!");
    }
    {
        Vec center(0.05, 0.95, 0.95);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.1);
        ngbsearch.set_boxsize(1.);

        Vec anchor(0.125, 0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == true,
                  "No hit on sphere intersects periodic copy case!");
    }

    return 0;
#else
    // No intersection
    {
        Vec center(0.5, 0.5);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.1);

        Vec anchor(0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == false,
                  "False hit on no intersection case!");
    }
    {
        Vec center(-0.5, -0.5);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.1);

        Vec anchor(0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == false,
                  "False hit on no intersection case!");
    }
    // Clear intersection
    {
        Vec center(0.25, 0.3);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.1);

        Vec anchor(0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == true,
                  "No hit on intersection case!");
    }
    {
        Vec center(0.25, -0.05);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.1);

        Vec anchor(0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == true,
                  "No hit on intersection case!");
    }
    // Sphere in box
    {
        Vec center(0.2, 0.2);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.01);

        Vec anchor(0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == true,
                  "No hit on sphere in box case!");
    }
    // Sphere touches box corner
    {
        Vec center(0.35, 0.35);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, sqrt(0.02));

        Vec anchor(0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == true,
                  "No hit on sphere touches box case!");
    }
    // Sphere intersects box side, but no corner
    {
        Vec center(0.125, 0.3);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.1);

        Vec anchor(0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == true,
                  "No hit on sphere intersects side case!");
    }
    {
        Vec center(-0.05, 0.125);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.1);

        Vec anchor(0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == true,
                  "No hit on sphere intersects side case!");
    }
    // Sphere intersects periodic copy
    {
        Vec center(0.95, 0.95);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.1);
        ngbsearch.set_boxsize(1.);

        Vec anchor(0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == true,
                  "No hit on sphere intersects periodic copy case!");
    }
    {
        Vec center(0.05, 0.95);
        StarParticle star(center);
        ClosestNgbSearch ngbsearch(&star, 0.1);
        ngbsearch.set_boxsize(1.);

        Vec anchor(0.125, 0.125);
        Box box(anchor, 0.25);
        TreeNode treenode(box);

        my_assert(ngbsearch.splitnode(&treenode) == true,
                  "No hit on sphere intersects periodic copy case!");
    }

    return 0;
#endif
}
