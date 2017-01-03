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
 * @file testYMLFile.cpp
 *
 * @brief Unit test for the YMLFile class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "YMLFile.hpp"
#include "myAssert.hpp"

/**
 * @brief Unit test for the YMLFile class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char** argv) {
    YMLFile file("test.yml");

    my_assert_simple(file.get_value<int>("test_integer1") == 42);
    my_assert_simple(file.get_value<int>("test_integer2") == 42);
    my_assert_simple(file.get_value<int>("test_integer3") == 42);
    my_assert_simple(file.get_value<double>("test_float") == 3.14);
    my_assert_simple(file.get_value<bool>("test_bool1") == true);
    my_assert_simple(file.get_value<bool>("test_bool2") == true);
    my_assert_simple(file.get_value<bool>("test_bool3") == true);
    my_assert_simple(file.get_value<bool>("test_bool4") == true);
    my_assert_simple(file.get_value<bool>("test_bool5") == false);
    my_assert_simple(file.get_value<bool>("test_bool6") == false);
    my_assert_simple(file.get_value<bool>("test_bool7") == false);
    my_assert_simple(file.get_value<bool>("test_bool8") == false);
    my_assert_simple(file.get_value<std::string>("test_string") ==
                     "This is a test string.");
    my_assert_simple(file.get_value<int>("test_group:test_group_member") == 42);
    my_assert_simple(file.get_value<std::string>(
                             "test_comments_group:test_comments_value") ==
                     "test comments string");

    my_assert_simple(
            file.get_value<int>(
                    "test_group2:test_group_group:test_group_group_member") ==
            42);

    // default values
    my_assert_simple(file.get_value<int>("not_in_file1", 42) == 42);
    my_assert_simple(file.get_value<double>("not_in_file2", 3.14) == 3.14);
    my_assert_simple(file.get_value<std::string>("not_in", "file?") == "file?");
    my_assert_simple(file.get_value<bool>("not_in_file3", true) == true);

    my_assert_simple(
            file.get_value<std::string>(
                    "test_group2:test_group_group:test_group_str_member",
                    "hello!") == "hello!");

    file.print_contents(std::cout);

    return 0;
}
