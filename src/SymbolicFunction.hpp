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
 * @file SymbolicFunction.hpp
 *
 * @brief Support for symbolic functions that are interpreted at runtime
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SYMBOLICFUNCIONS_HPP
#define SYMBOLICFUNCIONS_HPP

#define BOOST_SPIRIT_USE_PHOENIX_V3

#include <boost/config/warning_disable.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/qi.hpp>
#include <cmath>      // for acos, asin, atan, cbrt, cos, etc
#include <exception>  // for exception
#include <string>     // for basic_string, etc

namespace mathparser {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

/*! @brief Symbolic mathematical constants */
static struct constant_ : qi::symbols<char, double> {
    constant_() { this->add("pi", boost::math::constants::pi<double>()); }
} constant; /*!< @brief Symbolic mathematical constants */

/*! @brief Symbolic functions with a single argument */
static struct function_ : qi::symbols<char, double (*)(double)> {
    function_() {
        this->add("cos", static_cast<double (*)(double)>(&std::cos))(
                "sin", static_cast<double (*)(double)>(&std::sin))(
                "tan", static_cast<double (*)(double)>(&std::tan))(
                "sinh", static_cast<double (*)(double)>(&std::sinh))(
                "cosh", static_cast<double (*)(double)>(&std::cosh))(
                "tanh", static_cast<double (*)(double)>(&std::tanh))(
                "log", static_cast<double (*)(double)>(&std::log))(
                "exp", static_cast<double (*)(double)>(&std::exp))(
                "acos", static_cast<double (*)(double)>(&std::acos))(
                "asin", static_cast<double (*)(double)>(&std::asin))(
                "atan", static_cast<double (*)(double)>(&std::atan))(
                "log10", static_cast<double (*)(double)>(&std::log10))(
                "sqrt", static_cast<double (*)(double)>(&std::sqrt))(
                "cbrt", static_cast<double (*)(double)>(&cbrt));
    }
} function; /*!< @brief Symbolic functions with a single argument */

/*! @brief Symbolic power function */
struct power_ {

    /*! @brief Default result type */
    template <class> struct result;

    /*! @brief Phoenix 3 compliant result type */
    template <class F, typename X, typename Y> struct result<F(X, Y)> {
        /*! @brief The result is an X reference */
        typedef X& type;
    };

    /**
     * @brief Symbolic power operation
     *
     * @param x Base of the power
     * @param y Exponent of the power
     * @return Base to the power exponent
     */
    template <typename X, typename Y> X& operator()(X& x, Y y) const {
        x = std::pow(x, y);
        return x;
    }
};

/*! @brief Symbolic function wrapper around a function */
struct func_ {

    /*! @brief Default result type */
    template <class> struct result;

    /*! @brief Phoenix 3 compliant result type */
    template <class Q, typename F, typename X> struct result<Q(F, X)> {
        /*! @brief The result is an X reference */
        typedef X& type;
    };

    /**
     * @brief Symbolic function operation
     *
     * @param f Function to call
     * @param x Parameter passed on to the function
     * @return Result of the function call
     */
    template <typename F, typename X> X& operator()(F f, X& x) const {
        x = f(x);
        return x;
    }
};

/**
 * @brief Parser for symbolic mathematical expressions
 */
struct math : qi::grammar<std::string::const_iterator, double(),
                          ascii::space_type> {
    /**
     * @brief Constructor
     *
     * @param value_r Radius passed on to the expression
     * @param value_x x-coordinate passed on to the expression
     * @param value_y y-coordinate passed on to the expression
     * @param value_z z-coordinate passed on to the expression
     */
    math(double value_r = 1., double value_x = 1., double value_y = 1.,
         double value_z = 1.)
            : math::base_type(expr) {
        boost::phoenix::function<power_> power;
        boost::phoenix::function<func_> func;

        //        struct function_ function;
        //        struct constant_ constant;

        expr = term[qi::_val = qi::_1] >> *(('+' >> term[qi::_val += qi::_1]) |
                                            ('-' >> term[qi::_val -= qi::_1]));

        term = factor[qi::_val = qi::_1] >>
               *(('*' >> factor[qi::_val *= qi::_1]) |
                 ('/' >> factor[qi::_val /= qi::_1]));

        factor = arg[qi::_val = qi::_1] >>
                 *("**" >> factor[qi::_val = power(qi::_val, qi::_1)]);

        arg = qi::double_[qi::_val = qi::_1] | symbols[qi::_val = qi::_1] |
              ('(' >> expr[qi::_val = qi::_1] >> ')') |
              ('-' >> arg[qi::_val = -qi::_1]) | (constant[qi::_val = qi::_1]) |
              ((function >> '(' >> expr >>
                ')')[qi::_val = func(qi::_1, qi::_2)]);

        symbols = symbol_r[qi::_val = qi::_1] | symbol_x[qi::_val = qi::_1] |
                  symbol_y[qi::_val = qi::_1] | symbol_z[qi::_val = qi::_1];

        symbol_r = qi::lit('r')[qi::_val = value_r];
        symbol_x = qi::lit('x')[qi::_val = value_x];
        symbol_y = qi::lit('y')[qi::_val = value_y];
        symbol_z = qi::lit('z')[qi::_val = value_z];
    }

    /**
     * @brief Change the radius
     *
     * @param value New radius
     */
    void change_value_r(double value) {
        symbol_r = qi::lit('r')[qi::_val = value];
    }

    /**
     * @brief Change the x-coordinate
     *
     * @param value New x-coordinate
     */
    void change_value_x(double value) {
        symbol_x = qi::lit('x')[qi::_val = value];
    }

    /**
     * @brief Change the y-coordinate
     *
     * @param value New y-coordinate
     */
    void change_value_y(double value) {
        symbol_y = qi::lit('y')[qi::_val = value];
    }

    /**
     * @brief Change the z-coordinate
     *
     * @param value New z-coordinate
     */
    void change_value_z(double value) {
        symbol_z = qi::lit('z')[qi::_val = value];
    }

    /**@{*/
    /*! @brief Rules used to parse the expression */
    qi::rule<std::string::const_iterator, double(), ascii::space_type> expr,
            term, factor, arg, symbols, symbol_r, symbol_x, symbol_y, symbol_z;
    /**@}*/
};
}

/**
 * @brief Exception thrown when the string input to SymbolicFunction could not
 * be interpreted
 */
class badexpressionexception : public std::exception {
  public:
    /**
     * @brief Return a human-readable message
     *
     * @return A human-readable message
     */
    virtual const char* what() const throw() {
        return "A malformed expression was encountered";
    }
};

/**
 * @brief Mathematical function parsed from an input string
 *
 * Used to convert a string to a functional runtime function that can be used to
 * set up complex density/velocity/pressure... profiles without the need to
 * write hard code in SpecificICGenerator.
 */
class SymbolicFunction {
  private:
    /*! @brief Parser used to interpret strings and treat them as functions */
    mathparser::math _math;

    /*! @brief String that is being interpreted */
    std::string _str;

    /**
     * @brief Parse the given string into a functional mathparser
     *
     * @param str std::string to parse
     * @return Double precision floating point indicating success or failure (?)
     */
    double parse_expr(std::string str) {
        double result;
        std::string::const_iterator iter = str.begin();
        std::string::const_iterator end = str.end();
        bool r = phrase_parse(iter, end, _math, boost::spirit::ascii::space,
                              result);
        if(!r || iter != end) {
            result = nan("");
        }
        return result;
    }

  public:
    /**
     * @brief Constructor
     *
     * Create a SymbolicFunction based on the given std::string. Can be used as
     * if it  were an ordinary function (although it will be significantly
     * slower).
     *
     * @param str std::string representation of a mathematical function
     */
    SymbolicFunction(std::string str) {
        double result = parse_expr(str);
        if(result != result) {
            throw badexpressionexception();
        }
        _str = str;
    }

    /**
     * @brief Calculate the function value of the SymbolicFunction at the given
     * radius and/or coordinates
     *
     * @param r Radius at which to evaluate the function
     * @param x x-coordinate at which to evaluate the function
     * @param y y-coordinate at which to evaluate the function
     * @param z z-coordinate at which to evaluate the function
     * @return Value of the function at the given coordinate(s)
     */
    double operator()(double r = 1., double x = 1., double y = 1.,
                      double z = 1.) {
        _math.change_value_r(r);
        _math.change_value_x(x);
        _math.change_value_y(y);
        _math.change_value_z(z);
        return parse_expr(_str);
    }
};

#endif  // SYMBOLICFUNCIONS_HPP
