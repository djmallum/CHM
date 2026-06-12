//
// Created by Allum, Donovan on 2026-06-03.
//

#pragma once
#include <concepts>
#include "LinearAlgebra.hpp"
#include <string>

namespace SoilMoistureSolver::detail
{
class data;

template<class T,typename... Args>
concept ElementInterface = requires(const T& t,Args&&... args)
{
    { t->template make_module_data<data>(std::declval<std::string>(),std::forward<Args>(args)...)} -> std::same_as<data&>;
    { t->template get_module_data<data>(std::declval<std::string>())} -> std::same_as<data&>;
    { t->template edge_unit_normal<Vector_3>(std::declval<size_t>())} -> std::same_as<Vector_3>;
    { t->template edge_midpoint<Point_3>(std::declval<size_t>())} -> std::same_as<Point_3>;
    { t->downslope_dir() } -> std::same_as<Vector_3>;
    //{ t->soil_attribute(std::declval<std::string>())} -> std::floating_point;
    { t->neighbor(std::declval<size_t>())} -> std::same_as<T>;
    { t->cell_global_id } -> std::convertible_to<size_t>;
};

template<class T>
concept MeshInterface = requires(const T& t)
{
    { t->face(std::declval<int>()) } -> ElementInterface;
    { t->size_local_faces() } -> std::convertible_to<size_t>;
    { t->size_global_faces() } -> std::convertible_to<size_t>;
};

template<typename T>
concept FloatingPointUnderReference = std::floating_point<std::remove_cvref_t<T>>;

template<typename T>
concept Indexable = requires(T obj, size_t idx) {
    { obj[idx] } -> FloatingPointUnderReference;  // Return type convertible to size_t
};
}