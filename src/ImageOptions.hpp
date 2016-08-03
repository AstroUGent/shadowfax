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
 * @file ImageOptions.hpp
 *
 * @brief Options for images
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef IMAGEOPTIONS_HPP
#define IMAGEOPTIONS_HPP

// The ImageType is a combination of filetype info and colormap info (if
// applicable)
/**
  * @brief Image options
  */
enum ImageType {
    /*! @brief A binary grayscale image */
    BIN_PGM,
    /*! @brief An ascii grayscale image */
    ASC_PGM,
    /*! @brief A binary color image */
    BIN_PPM,
    /*! @brief An ascii color image */
    ASC_PPM,
    /*! @brief An image with the "jet" color scheme */
    CM_JET,
    /*! @brief An image with the "afmhot" color scheme */
    CM_AFMHOT = 8,
    /*! @brief An image with the "ocean" color scheme */
    CM_OCEAN = 12,
    /*! @brief An image with the "rdbu" color scheme */
    CM_RDBU = 16
};

#endif  // IMAGEOPTIONS_HPP
