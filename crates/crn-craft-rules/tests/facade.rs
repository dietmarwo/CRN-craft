//! The facade crate re-exports `crn_craft::rules`. Every public path of
//! crn-craft-rules 0.6.0 must still resolve, and resolve to the engine's own
//! item, so that user code and generated kernels and loops keep compiling.
//! The lists are the 0.6.0 rustdoc item index plus the root re-export
//! `Model` and the hidden `libm` re-export.

use std::any::TypeId;

fn id<T: ?Sized + 'static>() -> TypeId {
    TypeId::of::<T>()
}

fn id_of<T: 'static>(_: T) -> TypeId {
    TypeId::of::<T>()
}

macro_rules! same_types {
    ({ $($name:ident),* $(,)? }) => {
        $(assert_eq!(
            id::<crn_craft_rules::$name>(),
            id::<crn_craft::rules::$name>(),
            stringify!($name),
        );)*
    };
    ($module:ident::{ $($name:ident),* $(,)? }) => {
        $(assert_eq!(
            id::<crn_craft_rules::$module::$name>(),
            id::<crn_craft::rules::$module::$name>(),
            concat!(stringify!($module), "::", stringify!($name)),
        );)*
    };
}

macro_rules! same_functions {
    ({ $($name:ident),* $(,)? }) => {
        $(assert_eq!(
            id_of(crn_craft_rules::$name),
            id_of(crn_craft::rules::$name),
            stringify!($name),
        );)*
    };
    ($module:ident::{ $($name:ident),* $(,)? }) => {
        $(assert_eq!(
            id_of(crn_craft_rules::$module::$name),
            id_of(crn_craft::rules::$module::$name),
            concat!(stringify!($module), "::", stringify!($name)),
        );)*
    };
}

macro_rules! same_constants {
    ({ $($name:ident),* $(,)? }) => {
        $(assert_eq!(
            crn_craft_rules::$name,
            crn_craft::rules::$name,
            stringify!($name),
        );)*
    };
    ($module:ident::{ $($name:ident),* $(,)? }) => {
        $(assert_eq!(
            crn_craft_rules::$module::$name,
            crn_craft::rules::$module::$name,
            concat!(stringify!($module), "::", stringify!($name)),
        );)*
    };
}

struct Counter;

impl crn_craft_rules::observe::EventObserver for Counter {
    fn on_event(
        &self,
        _: &crn_craft_rules::observe::Event<'_>,
    ) -> crn_craft_rules::observe::ObserverAction {
        crn_craft_rules::observe::ObserverAction::Continue
    }
}

fn engine_observer<O: crn_craft::rules::observe::EventObserver>() -> bool {
    O::ENABLED
}

struct NoKernel;

impl crn_craft_rules::RulesKernel for NoKernel {
    fn fingerprint(&self) -> &str {
        ""
    }
    fn rates(&self, _: Option<usize>, _: &[u64], _: &[f64], _: &mut [f64]) -> Result<(), usize> {
        Ok(())
    }
    fn projection(&self, _: usize, _: &[u64]) -> f64 {
        0.0
    }
}

struct NoLoop;

impl crn_craft_rules::RulesLoop for NoLoop {
    fn fingerprint(&self) -> &str {
        ""
    }
    fn run(
        &self,
        _: &mut crn_craft_rules::LoopHost<'_, '_>,
        _: &mut [u64],
        _: &mut [f64],
    ) -> Result<Option<f64>, crn_craft_rules::LoopHalt> {
        Ok(None)
    }
}

type Outcome = Result<crn_craft_rules::RulesResult, crn_craft::diagnostic::Diagnostic>;
type KernelEntry = fn(&crn_craft_rules::Prepared, &crn_craft_rules::Control, &NoKernel) -> Outcome;
type LoopEntry = fn(&crn_craft_rules::Prepared, &crn_craft_rules::Control, &NoLoop) -> Outcome;

#[test]
fn root_items_resolve_to_the_engine() {
    same_types!({
        Budget,
        CheckVersion,
        Control,
        FinalInstances,
        InstanceRecord,
        LoopHalt,
        Method,
        Model,
        Outputs,
        PathSeries,
        Prepared,
        Record,
        ResultVersion,
        RulesCheck,
        RulesResult,
        RulesWaitingTime,
        RunRequest,
        RunVersion,
        Statistics,
        Termination,
        Time,
        TypeHeader,
    });
    assert_eq!(
        id::<crn_craft_rules::LoopHost<'static, 'static>>(),
        id::<crn_craft::rules::LoopHost<'static, 'static>>(),
    );
    assert_eq!(
        id::<dyn crn_craft_rules::RulesKernel>(),
        id::<dyn crn_craft::rules::RulesKernel>(),
    );
    assert_eq!(
        id::<dyn crn_craft_rules::RulesLoop>(),
        id::<dyn crn_craft::rules::RulesLoop>(),
    );
    same_functions!({
        check_model,
        check_model_file,
        describe,
        execute,
        generate_rules_kernel,
        generate_rules_loop,
        prepare,
        schema,
        solve,
        solve_files,
        stream,
    });
    // A kernel or loop written against the facade traits runs on the engine.
    let _: [KernelEntry; 2] = [
        crn_craft_rules::execute_with_rules_kernel,
        crn_craft::rules::execute_with_rules_kernel,
    ];
    let _: [LoopEntry; 2] = [
        crn_craft_rules::execute_with_rules_loop,
        crn_craft::rules::execute_with_rules_loop,
    ];
    assert_eq!(
        id_of(crn_craft_rules::execute_with_observer::<Counter>),
        id_of(crn_craft::rules::execute_with_observer::<Counter>),
    );
    same_constants!({
        MAX_REPLICATES,
        MAX_REPORT_TIMES,
        MAX_RESULT_BYTES,
        MAX_WORKERS,
        MAX_WORKSPACE_BYTES,
        SAMPLER,
    });
}

#[test]
fn public_modules_resolve_to_the_engine() {
    same_types!(atp::{ AtpSpec });
    same_types!(continuous::{ ContinuousSpec, OnBound });
    same_constants!(continuous::{ MAX_VARIABLES });
    same_types!(model::{
        Initial,
        Model,
        ModelVersion,
        Observable,
        Pattern,
        Product,
        Reduce,
        Rule,
        SpeciesType,
        Value,
    });
    same_constants!(model::{
        MAX_ATTRIBUTES,
        MAX_COUNT,
        MAX_EXPRESSION_BYTES,
        MAX_INITIAL,
        MAX_NAME_BYTES,
        MAX_OBSERVABLES,
        MAX_PARAMETERS,
        MAX_PRODUCTS,
        MAX_RULES,
        MAX_TYPES,
    });
    same_types!(observe::{
        Episode,
        LogEntry,
        LogSpec,
        Observation,
        ObserverAction,
        TrackerResult,
        TrackerSpec,
    });
    assert_eq!(
        id::<crn_craft_rules::observe::Event<'static>>(),
        id::<crn_craft::rules::observe::Event<'static>>(),
    );
    same_functions!(observe::{ log_entry });
    // An observer written against the facade trait is an engine observer.
    assert!(engine_observer::<Counter>());
}

#[test]
fn hidden_libm_reexport_serves_generated_code() {
    // Generated stop predicates and rates call `crn_craft_rules::libm::*`.
    assert_eq!(crn_craft_rules::libm::pow(2.0, 10.0), 1024.0);
    assert_eq!(
        id_of(crn_craft_rules::libm::pow),
        id_of(crn_craft::rules::libm::pow)
    );
    for f in [
        crn_craft_rules::libm::ceil,
        crn_craft_rules::libm::exp,
        crn_craft_rules::libm::floor,
        crn_craft_rules::libm::log,
        crn_craft_rules::libm::sqrt,
    ] {
        assert!(f(1.0).is_finite());
    }
}
