//
// Created by Allum, Donovan on 2026-06-03.
//

#pragma once
#include <concepts>
#include "LinearAlgebra.hpp"
#include <string>

namespace SoilMoistureSolver::detail
{
class Data;
template<class T>
concept Element = requires(const T& t)
{
    { t->template make_module_data<Data>(std::declval<std::string>())} -> std::same_as<Data&>;
    { t->template get_module_data<Data>(std::declval<std::string>())} -> std::same_as<Data&>;
    { t->template edge_unit_normal<Vector_3>(std::declval<size_t>())} -> std::same_as<Vector_3>;
    { t->template edge_midpoint<Point_3>(std::declval<size_t>())} -> std::same_as<Point_3>;
    { t->soil_attribute(std::declval<std::string>())} -> std::floating_point;
};

template<class T>
concept Mesh = requires(const T& t)
{
    { t->face(std::declval<int>()) } -> Element;
    { t->size_local_faces() } -> std::convertible_to<size_t>;
    { t->size_global_faces() } -> std::convertible_to<size_t>;
};

template<typename T>
concept Indexable = requires(T obj, size_t idx) {
    { obj[idx] } -> std::floating_point;  // Return type convertible to size_t
};
}