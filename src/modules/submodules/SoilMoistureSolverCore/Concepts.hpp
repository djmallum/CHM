//
// Created by Allum, Donovan on 2026-06-03.
//

#pragma once
#include <concepts>
#include "LinearAlgebra.hpp"
#include <string>

namespace SoilMoistureSolver::detail
{
template<class T,class Data,typename... Args>
concept ElementInterface = requires(const T& t,Args&&... args,const std::string& string_ref)
{
    { t->template make_module_data<Data>(string_ref,std::forward<Args>(args)...)} -> std::same_as<Data&>;
    { t->template get_module_data<Data>(string_ref,std::declval<std::string>())} -> std::same_as<Data&>;
    { t->template edge_unit_normal<Vector_3>(std::declval<size_t>())} -> std::same_as<Vector_3>;
    { t->template edge_midpoint<Point_3>(std::declval<size_t>())} -> std::same_as<Point_3>;
    { t->downslope_dir() } -> std::same_as<Vector_3>;
    //{ t->soil_attribute(std::declval<std::string>())} -> std::floating_point;
    { t->neighbor(std::declval<size_t>())} -> std::same_as<T>;
    { t->cell_global_id } -> std::convertible_to<size_t>;
};

template<class T, class Data>
concept MeshInterface = requires(const T& t)
{
    { t->face(std::declval<int>()) } -> ElementInterface<Data>;
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