#pragma once

// Standard
#include <concepts>
#include <ranges>
#include <string>
#include <tuple>
#include <type_traits>
#include <variant>
#include <vector>

// Internal
#include "EvaluationContext.hpp"

namespace detail::meta {
template <typename T, typename VARIANT_T>
struct is_variant_member;

template <typename T, typename... ALL_T>
struct is_variant_member<T, std::variant<ALL_T...>>
    : public std::disjunction<std::is_same<T, ALL_T>...> {};

template <typename T, typename VARIANT>
static constexpr bool is_variant_member_v = is_variant_member<T, VARIANT>::value;

// a simple list of types
template <typename... Ts>
struct list {};

// append two lists
template <typename A, typename B>
struct append;
template <typename... As, typename... Bs>
struct append<list<As...>, list<Bs...>> {
    using type = list<As..., Bs...>;
};
template <typename A, typename B>
using append_t = typename append<A, B>::type;

// concatenate any number of lists
template <typename... Lists>
struct concat;
template <>
struct concat<> {
    using type = list<>;
};
template <typename L>
struct concat<L> {
    using type = L;
};
template <typename L1, typename L2, typename... Rest>
struct concat<L1, L2, Rest...> : concat<append_t<L1, L2>, Rest...> {};
template <typename... Ls>
using concat_t = typename concat<Ls...>::type;

// given L and list<Rs...>, build list< list<L,Rs>... >
template <typename L, typename RList>
struct cross_pair;
template <typename L, typename... Rs>
struct cross_pair<L, list<Rs...>> {
    using type = list<list<L, Rs>...>;
};
template <typename L, typename RList>
using cross_pair_t = typename cross_pair<L, RList>::type;

// cross-product of two lists
template <typename LList, typename RList>
struct cross;
template <typename... Ls, typename... Rs>
struct cross<list<Ls...>, list<Rs...>> {
    using type = concat_t<cross_pair_t<Ls, list<Rs...>>...>;
};
template <typename LList, typename RList>
using cross_t = typename cross<LList, RList>::type;

// map an MPL list of pairs< G,EL > → a list of EvaluationContext<G,EL...>
template <typename Pair>
struct to_context;
template <typename G, typename... Es>
struct to_context<list<G, list<Es...>>> {
    using type = pipelines::detect::EvaluationContext<G, Es...>;
};

template <typename List>
struct variant_of;
template <typename... Ts>
struct variant_of<list<Ts...>> {
    using type = std::variant<Ts...>;
};
template <typename L>
using variant_of_t = typename variant_of<L>::type;

template <typename P>
using to_context_t = typename to_context<P>::type;

template <typename PairList>
struct transform_contexts;
template <typename... Ps>
struct transform_contexts<list<Ps...>> {
    using type = list<to_context_t<Ps>...>;
};
template <typename PL>
using transform_contexts_t = typename transform_contexts<PL>::type;

template <typename GList, typename EList>
using contexts_list_t = transform_contexts_t<cross_t<GList, EList>>;

template <typename GList, typename EList>
using contexts_variant_t = variant_of_t<contexts_list_t<GList, EList>>;

// primary template for prefixes:
template <typename List>
struct prefixes;

// base case: an empty list has exactly one prefix, the empty list
template <>
struct prefixes<list<>> {
    using type = list<list<>>;
};

// inductive case: peel off First, compute prefixes of the rest,
// then “prepend” First to each of those, and finally stick on the empty prefix
template <typename First, typename... Rest>
struct prefixes<list<First, Rest...>> {
   private:
    // prefixes of the tail
    using tail = typename prefixes<list<Rest...>>::type;

    // helper to prepend First to a single prefix
    template <typename P>
    struct prepend;

    template <typename... Ps>
    struct prepend<list<Ps...>> {
        using type = list<First, Ps...>;
    };

    // map prepend<> over every prefix in tail
    template <typename ListOf>
    struct map_prepend;

    template <typename... Ps>
    struct map_prepend<list<Ps...>> {
        using type = list<typename prepend<Ps>::type...>;
    };

    using extended = typename map_prepend<tail>::type;

   public:
    // result = empty-prefix + all the extended ones
    using type = typename concat<list<list<>>, extended>::type;
};

// helper alias:
template <typename L>
using prefixes_t = typename prefixes<L>::type;

template <typename InitialVariant, typename... Steps>
struct make_step_variants;

//------------------------------------------------------------------------------
// Flatten nested variants: type_list_from_variant
//------------------------------------------------------------------------------

// Fallback: non-variant type becomes single-element list

template <typename T>
struct type_list_from_variant {
    using type = list<T>;
};

template <typename T>
using type_list_from_variant_t = typename type_list_from_variant<T>::type;

// If T is std::variant<Us...>, extract Us... into list

template <typename... Us>
struct type_list_from_variant<std::variant<Us...>> {
    using type = concat_t<type_list_from_variant_t<Us>...>;
};

template <typename R>
    requires(!std::same_as<R, std::string>) && std::ranges::range<R>
struct type_list_from_variant<R> {
    // peel off the container, *then* flatten any variants in the element-type
    using type = type_list_from_variant_t<std::ranges::range_value_t<R>>;
};

//------------------------------------------------------------------------------
// Step application: apply_t
//------------------------------------------------------------------------------

template <class Step, class Ctx, class = void>
struct apply_impl {
    using type = Ctx;
};

template <class Step, class Ctx>
struct apply_impl<Step, Ctx, std::void_t<std::invoke_result_t<Step, Ctx>>> {
    using type = std::invoke_result_t<Step, Ctx>;
};

template <class Step, class Ctx>
using apply_t = typename apply_impl<Step, Ctx>::type;

//------------------------------------------------------------------------------
// transform_variant: apply Step over each Variant alternative and flatten
//------------------------------------------------------------------------------

template <class Step, class Variant>
struct transform_variant;

// template <class Step, typename... Ctxs>
// struct transform_variant<Step, std::variant<Ctxs...>> {
//     using flat_list = concat_t<type_list_from_variant_t<apply_t<Step, Ctxs>>...>;

//     using type = variant_of_t<flat_list>;
// };
//
template <class Step, typename... Ctxs>
struct transform_variant<Step, std::variant<Ctxs...>> {
   private:
    // primary: for Ctx where Step(ctx) is ill-formed
    template <typename Ctx, bool Inv = std::is_invocable_v<Step, Ctx>>
    struct step_list_impl {
        using type = list<Ctx>;  // carry it through unchanged
    };

    // specialization: only instantiated when Invocable == true
    template <typename Ctx>
    struct step_list_impl<Ctx, true> {
        using type = type_list_from_variant_t<std::invoke_result_t<Step, Ctx>>;
    };

    template <typename Ctx>
    using step_list_t = typename step_list_impl<Ctx>::type;

    using flat_list = concat_t<step_list_t<Ctxs>...>;

   public:
    using type = variant_of_t<flat_list>;
};

template <class Step, class Variant>
using transform_variant_t = typename transform_variant<Step, Variant>::type;

//------------------------------------------------------------------------------
// apply_all_steps: fold over Steps, keeping variant flattened
//------------------------------------------------------------------------------

template <typename Variant, typename... Steps>
struct apply_all_steps;

template <typename Variant>
struct apply_all_steps<Variant> {
    using type = Variant;
};

template <typename Variant, typename Step, typename... Rest>
struct apply_all_steps<Variant, Step, Rest...> {
    using next = transform_variant_t<Step, Variant>;
    using type = typename apply_all_steps<next, Rest...>::type;
};

template <typename Variant, typename... Steps>
using apply_all_steps_t = typename apply_all_steps<Variant, Steps...>::type;

//------------------------------------------------------------------------------
// make_step_variants: tuple of successive variant stages (flattened)
//------------------------------------------------------------------------------

template <typename InitialVariant, typename... Steps>
struct make_step_variants;

template <typename InitialVariant>
struct make_step_variants<InitialVariant> {
    using type = std::tuple<InitialVariant>;
};

template <typename InitialVariant, typename Step, typename... Rest>
struct make_step_variants<InitialVariant, Step, Rest...> {
    using NextVar = transform_variant_t<Step, InitialVariant>;
    using Tail = typename make_step_variants<NextVar, Rest...>::type;
    using type =
        decltype(std::tuple_cat(std::declval<std::tuple<InitialVariant>>(), std::declval<Tail>()));
};

template <typename InitialVariant, typename... Steps>
using make_step_variants_t = typename make_step_variants<InitialVariant, Steps...>::type;
template <typename Tuple>
struct tuple_of_vector_from_tuple;

template <typename... Ts>
struct tuple_of_vector_from_tuple<std::tuple<Ts...>> {
    using type = std::tuple<std::vector<Ts>...>;
};

template <typename Tuple>
using tuple_of_vector_from_tuple_t = typename tuple_of_vector_from_tuple<Tuple>::type;

}  // namespace detail::meta
