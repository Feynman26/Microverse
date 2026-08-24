extends RefCounted
class_name ExperimentAnalytics

# Analytics are deliberately post hoc. They accept immutable run dictionaries
# and return hypotheses/summary statistics; they never alter physiology, RNG,
# environment, heredity or survival.

static func summarize_batch(batch_or_runs) -> Dictionary:
	var runs: Array = batch_or_runs.get("runs", []) if batch_or_runs is Dictionary else batch_or_runs
	assert(not runs.is_empty())
	var final_populations: Array = []
	var max_populations: Array = []
	var mutation_counts: Array = []
	var division_counts: Array = []
	var termination_counts: Dictionary = {}
	var extinction_count: int = 0
	var seeds: Array = []
	for run_variant in runs:
		var run: Dictionary = run_variant
		seeds.append(int(run.get("seed", 0)))
		final_populations.append(float(run.get("final_population", 0)))
		max_populations.append(float(run.get("max_population", 0)))
		mutation_counts.append(float(run.get("mutation_events", 0)))
		division_counts.append(float(run.get("division_events", 0)))
		var termination: String = String(run.get("termination_reason", "unknown"))
		termination_counts[termination] = int(termination_counts.get(termination, 0)) + 1
		if int(run.get("final_population", 0)) == 0:
			extinction_count += 1
	return {
		"replicate_count": runs.size(),
		"seeds": seeds,
		"extinction_count": extinction_count,
		"extinction_frequency": float(extinction_count) / float(runs.size()),
		"termination_counts": termination_counts,
		"final_population": _distribution(final_populations),
		"max_population": _distribution(max_populations),
		"mutation_events": _distribution(mutation_counts),
		"division_events": _distribution(division_counts)
	}

static func detect_candidates(run: Dictionary) -> Dictionary:
	var trajectory: Array = run.get("trajectory", [])
	assert(not trajectory.is_empty())
	var populations: Array = []
	for sample_variant in trajectory:
		populations.append(int(sample_variant.get("population", 0)))
	var initial_population: int = int(populations[0])
	var final_population: int = int(populations[-1])
	var max_population: int = 0
	var min_population: int = 2147483647
	for population_variant in populations:
		var population: int = int(population_variant)
		max_population = maxi(max_population, population)
		min_population = mini(min_population, population)

	var persistent_polymorphism: bool = _persistent_polymorphism(trajectory)
	var turning_points: int = _turning_points(populations)
	var expansion: bool = initial_population > 0 and final_population >= initial_population * 2
	var bottleneck: bool = max_population >= 4 and min_population <= maxi(1, int(floor(float(max_population) * 0.25)))
	var extinct: bool = final_population == 0
	var rescue: bool = _ecological_rescue(populations)
	var mutation_shift: Dictionary = _mutation_rate_shift(trajectory)
	var genome_shift: Dictionary = _genome_size_shift(trajectory)
	var resource_novelty: Dictionary = _closed_resource_novelty(run)

	return {
		"rapid_persistent_lineage_expansion": {
			"candidate": expansion and _tail_noncollapsing(populations),
			"initial_population": initial_population,
			"final_population": final_population
		},
		"persistent_polymorphism": {
			"candidate": persistent_polymorphism,
			"evidence_samples": mini(3, trajectory.size())
		},
		"population_bottleneck": {
			"candidate": bottleneck,
			"minimum_population": min_population,
			"maximum_population": max_population
		},
		"complete_extinction": {
			"candidate": extinct,
			"termination_reason": String(run.get("termination_reason", "unknown"))
		},
		"new_resource_utilization": resource_novelty,
		"loss_of_former_resource_or_pathway": {
			"candidate": false,
			"evidence_available": false,
			"reason": "requires sampled reaction/transport usage history"
		},
		"new_secretion_uptake_coupling": {
			"candidate": false,
			"evidence_available": false,
			"reason": "requires sampled secretion/transport coupling history"
		},
		"genome_expansion_or_reduction": genome_shift,
		"mutation_rate_shift": mutation_shift,
		"recurring_population_cycles": {
			"candidate": turning_points >= 3,
			"turning_points": turning_points
		},
		"ecological_rescue_after_shock": {
			"candidate": rescue,
			"minimum_population": min_population,
			"final_population": final_population
		}
	}

static func compare_paired_fork(fork_result: Dictionary, metric: String = "final_population") -> Dictionary:
	var arms: Dictionary = fork_result.get("arms", {})
	assert(arms.size() >= 2)
	var values: Dictionary = {}
	var names: Array = arms.keys()
	names.sort()
	for name_variant in names:
		var name: String = String(name_variant)
		values[name] = float(arms[name_variant].get(metric, 0.0))
	var first_name: String = String(names[0])
	var second_name: String = String(names[1])
	return {
		"metric": metric,
		"values": values,
		"first_minus_second": float(values[first_name]) - float(values[second_name]),
		"fork_checksum": float(fork_result.get("fork_checksum", 0.0)),
		"prefix_checksums": fork_result.get("prefix_checksums", {}).duplicate(true)
	}

static func _persistent_polymorphism(trajectory: Array) -> bool:
	if trajectory.size() < 2:
		return false
	var start: int = maxi(0, trajectory.size() - 3)
	for index in range(start, trajectory.size()):
		if int(trajectory[index].get("genotype_count", 0)) < 2:
			return false
	return true

static func _tail_noncollapsing(populations: Array) -> bool:
	if populations.size() < 2:
		return false
	var start: int = maxi(0, populations.size() - 3)
	var tail_min: int = 2147483647
	for index in range(start, populations.size()):
		tail_min = mini(tail_min, int(populations[index]))
	return tail_min >= maxi(1, int(populations[0]))

static func _turning_points(populations: Array) -> int:
	if populations.size() < 3:
		return 0
	var result: int = 0
	var previous_direction: int = 0
	for index in range(1, populations.size()):
		var delta: int = int(populations[index]) - int(populations[index - 1])
		var direction: int = 1 if delta > 0 else (-1 if delta < 0 else 0)
		if direction == 0:
			continue
		if previous_direction != 0 and direction != previous_direction:
			result += 1
		previous_direction = direction
	return result

static func _ecological_rescue(populations: Array) -> bool:
	if populations.size() < 3:
		return false
	var initial: int = int(populations[0])
	var minimum: int = 2147483647
	var minimum_index: int = 0
	for index in range(1, populations.size() - 1):
		if int(populations[index]) < minimum:
			minimum = int(populations[index])
			minimum_index = index
	if minimum_index <= 0:
		return false
	var final_population: int = int(populations[-1])
	return minimum < initial and final_population > minimum and final_population >= initial

# Prefer M10's mechanistic expected per-copy error distribution when available.
# For archived pre-M10 run dictionaries, fall back to the older event-count
# heuristic rather than rewriting historical evidence.
static func _mutation_rate_shift(trajectory: Array) -> Dictionary:
	if trajectory.size() < 2:
		return {"candidate": false, "evidence_available": false, "reason": "insufficient samples"}
	var first_dynamics: Dictionary = trajectory[0].get("mutation_dynamics", {})
	var last_dynamics: Dictionary = trajectory[-1].get("mutation_dynamics", {})
	if not first_dynamics.is_empty() and not last_dynamics.is_empty():
		var first_dist: Dictionary = first_dynamics.get("point_error_rate_per_gene", {})
		var last_dist: Dictionary = last_dynamics.get("point_error_rate_per_gene", {})
		if int(first_dist.get("count", 0)) > 0 and int(last_dist.get("count", 0)) > 0:
			var initial_mean: float = float(first_dist.get("mean", 0.0))
			var final_mean: float = float(last_dist.get("mean", 0.0))
			var delta: float = final_mean - initial_mean
			return {
				"candidate": absf(delta) > 1e-12,
				"evidence_available": true,
				"initial_mean_point_error_rate_per_gene": initial_mean,
				"final_mean_point_error_rate_per_gene": final_mean,
				"change": delta,
				"source": "replication-derived expected copy error"
			}
	if trajectory.size() < 4:
		return {"candidate": false, "evidence_available": false, "reason": "pre-M10 trajectory has insufficient samples"}
	var midpoint: int = trajectory.size() / 2
	var first_start: int = int(trajectory[0].get("mutation_events", 0))
	var first_end: int = int(trajectory[midpoint - 1].get("mutation_events", 0))
	var second_start: int = int(trajectory[midpoint].get("mutation_events", 0))
	var second_end: int = int(trajectory[-1].get("mutation_events", 0))
	var first_delta: int = maxi(0, first_end - first_start)
	var second_delta: int = maxi(0, second_end - second_start)
	return {
		"candidate": second_delta >= maxi(2, first_delta * 2),
		"evidence_available": true,
		"source": "legacy realized-event heuristic"
	}

static func _genome_size_shift(trajectory: Array) -> Dictionary:
	if trajectory.size() < 2:
		return {"candidate": false, "evidence_available": false, "reason": "insufficient samples"}
	var first_dynamics: Dictionary = trajectory[0].get("mutation_dynamics", {})
	var last_dynamics: Dictionary = trajectory[-1].get("mutation_dynamics", {})
	if first_dynamics.is_empty() or last_dynamics.is_empty():
		return {
			"candidate": false,
			"evidence_available": false,
			"reason": "trajectory predates sampled M10 genome-size analytics"
		}
	var first_genes: Dictionary = first_dynamics.get("gene_count", {})
	var last_genes: Dictionary = last_dynamics.get("gene_count", {})
	var first_units: Dictionary = first_dynamics.get("replication_units", {})
	var last_units: Dictionary = last_dynamics.get("replication_units", {})
	if int(first_genes.get("count", 0)) <= 0 or int(last_genes.get("count", 0)) <= 0:
		return {"candidate": false, "evidence_available": false, "reason": "no living genomes in comparison samples"}
	var initial_mean_genes: float = float(first_genes.get("mean", 0.0))
	var final_mean_genes: float = float(last_genes.get("mean", 0.0))
	var initial_units: float = float(first_units.get("mean", 0.0))
	var final_units: float = float(last_units.get("mean", 0.0))
	return {
		"candidate": absf(final_units - initial_units) > 1e-12,
		"evidence_available": true,
		"initial_mean_gene_count": initial_mean_genes,
		"final_mean_gene_count": final_mean_genes,
		"gene_count_change": final_mean_genes - initial_mean_genes,
		"initial_mean_replication_units": initial_units,
		"final_mean_replication_units": final_units,
		"replication_unit_change": final_units - initial_units
	}

static func _closed_resource_novelty(run: Dictionary) -> Dictionary:
	var environment: Dictionary = run.get("environment", {})
	if String(environment.get("mode", "closed")) != "closed":
		return {
			"candidate": false,
			"evidence_available": false,
			"reason": "resource-total consumption is confounded by active boundary replacement"
		}
	var trajectory: Array = run.get("trajectory", [])
	if trajectory.size() < 2:
		return {"candidate": false, "evidence_available": false, "reason": "insufficient samples"}
	var first_resources: Dictionary = trajectory[0].get("resources", {})
	var last_resources: Dictionary = trajectory[-1].get("resources", {})
	var consumed: Dictionary = {}
	for field_variant in first_resources.keys():
		var field_name: String = String(field_variant)
		var delta: float = float(first_resources[field_variant]) - float(last_resources.get(field_variant, 0.0))
		if delta > 1e-9:
			consumed[field_name] = delta
	return {
		"candidate": not consumed.is_empty(),
		"evidence_available": true,
		"consumed_fields": consumed
	}

static func _distribution(values: Array) -> Dictionary:
	assert(not values.is_empty())
	var sorted: Array = values.duplicate()
	sorted.sort()
	return {
		"min": float(sorted[0]),
		"q25": _quantile_sorted(sorted, 0.25),
		"median": _quantile_sorted(sorted, 0.50),
		"q75": _quantile_sorted(sorted, 0.75),
		"max": float(sorted[-1]),
		"mean": _mean(sorted)
	}

static func _quantile_sorted(sorted: Array, fraction: float) -> float:
	assert(not sorted.is_empty() and fraction >= 0.0 and fraction <= 1.0)
	if sorted.size() == 1:
		return float(sorted[0])
	var position: float = fraction * float(sorted.size() - 1)
	var low: int = int(floor(position))
	var high: int = int(ceil(position))
	if low == high:
		return float(sorted[low])
	var weight: float = position - float(low)
	return lerpf(float(sorted[low]), float(sorted[high]), weight)

static func _mean(values: Array) -> float:
	var total: float = 0.0
	for value_variant in values:
		total += float(value_variant)
	return total / float(values.size())
