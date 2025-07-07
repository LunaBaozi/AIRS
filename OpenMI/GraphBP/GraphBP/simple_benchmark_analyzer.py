#!/usr/bin/env python3
"""
Simple Snakemake Benchmark Analyzer
Lightweight version that doesn't require external plotting libraries
"""

import os
import json
import glob
import csv
from datetime import datetime
from pathlib import Path

def parse_benchmark_file(filepath):
    """Parse individual Snakemake benchmark file"""
    try:
        data = {}
        with open(filepath, 'r') as f:
            lines = f.readlines()
            
        if len(lines) < 2:
            return None
            
        # Parse header and data
        headers = lines[0].strip().split('\t')
        values = lines[1].strip().split('\t')
        
        # Create data dictionary
        for i, header in enumerate(headers):
            if i < len(values):
                try:
                    data[header] = float(values[i])
                except ValueError:
                    data[header] = values[i]
        
        return {
            'rule_name': filepath.stem,
            'runtime_seconds': data.get('s', 0),
            'runtime_minutes': data.get('s', 0) / 60,
            'max_rss_mb': data.get('max_rss', 0) / 1024,  # Convert KB to MB
            'max_vms_mb': data.get('max_vms', 0) / 1024,
            'cpu_time': data.get('cpu_time', 0),
            'io_in': data.get('io_in', 0),
            'io_out': data.get('io_out', 0),
            'mean_load': data.get('mean_load', 0)
        }
    except Exception as e:
        print(f"Warning: Could not parse {filepath}: {e}")
        return None

def analyze_benchmarks(benchmark_dir="benchmarks", experiment=None):
    """Analyze benchmark files and generate report"""
    benchmark_dir = Path(benchmark_dir)
    
    # If experiment is specified, look in experiment subfolder
    if experiment:
        benchmark_dir = benchmark_dir / f"experiment_{experiment}"
        print(f"Analyzing experiment {experiment} benchmarks from {benchmark_dir}")
    
    if not benchmark_dir.exists():
        print(f"Benchmark directory '{benchmark_dir}' not found")
        return
    
    # Find all benchmark files
    benchmark_files = list(benchmark_dir.glob("*.txt"))
    
    if not benchmark_files:
        print(f"No benchmark files found in {benchmark_dir}")
        print("Make sure your Snakefile has 'benchmark:' directives in the rules")
        return
    
    print(f"Found {len(benchmark_files)} benchmark files")
    
    # Parse all benchmark data
    benchmark_data = []
    for file in benchmark_files:
        parsed = parse_benchmark_file(file)
        if parsed:
            benchmark_data.append(parsed)
    
    if not benchmark_data:
        print("No valid benchmark data found")
        return
    
    # Sort by runtime
    benchmark_data.sort(key=lambda x: x['runtime_seconds'], reverse=True)
    
    # Calculate totals
    total_runtime = sum(item['runtime_seconds'] for item in benchmark_data)
    total_memory = sum(item['max_rss_mb'] for item in benchmark_data)
    
    # Generate report
    print("\n" + "="*60)
    print("SNAKEMAKE PIPELINE BENCHMARK ANALYSIS")
    print("="*60)
    print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Rules analyzed: {len(benchmark_data)}")
    print(f"Total runtime: {total_runtime:.2f}s ({total_runtime/60:.2f} min)")
    print(f"Total memory: {total_memory:.2f} MB ({total_memory/1024:.2f} GB)")
    
    print("\n" + "-"*60)
    print("PERFORMANCE BY RULE:")
    print("-"*60)
    print(f"{'Rule Name':<25} {'Time (min)':<12} {'Memory (MB)':<12} {'% of Total':<10}")
    print("-"*60)
    
    for item in benchmark_data:
        percentage = (item['runtime_seconds'] / total_runtime) * 100
        print(f"{item['rule_name']:<25} {item['runtime_minutes']:<12.2f} {item['max_rss_mb']:<12.2f} {percentage:<10.1f}%")
    
    # Performance insights
    slowest_rule = benchmark_data[0]
    memory_hungry = max(benchmark_data, key=lambda x: x['max_rss_mb'])
    
    print("\n" + "-"*60)
    print("PERFORMANCE INSIGHTS:")
    print("-"*60)
    print(f"Slowest rule: {slowest_rule['rule_name']} ({slowest_rule['runtime_minutes']:.2f} min)")
    print(f"Most memory-intensive: {memory_hungry['rule_name']} ({memory_hungry['max_rss_mb']:.2f} MB)")
    
    # Find bottlenecks
    bottlenecks = [item for item in benchmark_data if item['runtime_seconds'] > total_runtime * 0.1]
    if bottlenecks:
        print(f"Rules using >10% of total time: {len(bottlenecks)}")
        for bottleneck in bottlenecks:
            percentage = (bottleneck['runtime_seconds'] / total_runtime) * 100
            print(f"   - {bottleneck['rule_name']}: {percentage:.1f}%")
    
    # CPU efficiency analysis
    cpu_efficient_rules = []
    for item in benchmark_data:
        if item['cpu_time'] > 0 and item['runtime_seconds'] > 0:
            efficiency = (item['cpu_time'] / item['runtime_seconds']) * 100
            item['cpu_efficiency'] = efficiency
            if efficiency < 50:
                cpu_efficient_rules.append(item)
    
    if cpu_efficient_rules:
        print(f"Rules with low CPU utilization (<50%):")
        for rule in cpu_efficient_rules:
            print(f"   - {rule['rule_name']}: {rule['cpu_efficiency']:.1f}%")
    
    print("\n" + "-"*60)
    print("OPTIMIZATION RECOMMENDATIONS:")
    print("-"*60)
    
    if slowest_rule['runtime_minutes'] > 5:
        print(f"• Focus optimization efforts on '{slowest_rule['rule_name']}' (biggest time sink)")
    
    if memory_hungry['max_rss_mb'] > 2000:  # > 2GB
        print(f"• Monitor memory usage for '{memory_hungry['rule_name']}' ({memory_hungry['max_rss_mb']:.0f} MB)")
    
    if len(bottlenecks) > 1:
        print(f"• Consider parallelizing or optimizing the {len(bottlenecks)} slowest rules")
    
    if cpu_efficient_rules:
        print(f"• Investigate I/O bottlenecks in rules with low CPU utilization")
    
    # Save CSV summary
    output_dir = Path("benchmark_results")
    if experiment:
        output_dir = output_dir / f"experiment_{experiment}"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    csv_file = output_dir / "benchmark_summary.csv"
    with open(csv_file, 'w', newline='') as f:
        if benchmark_data:
            writer = csv.DictWriter(f, fieldnames=benchmark_data[0].keys())
            writer.writeheader()
            writer.writerows(benchmark_data)
    
    print(f"\nDetailed data saved to: {csv_file}")
    print("="*60)

def analyze_stats_file(stats_file="benchmarks/execution_stats.json"):
    """Analyze Snakemake execution statistics"""
    if not os.path.exists(stats_file):
        print(f"Stats file not found: {stats_file}")
        print("   (This is normal for newer Snakemake versions)")
        return
    
    try:
        with open(stats_file, 'r') as f:
            stats = json.load(f)
        
        print(f"\nSNAKEMAKE EXECUTION STATISTICS:")
        print("-"*40)
        
        if 'total_runtime' in stats:
            print(f"Total pipeline runtime: {stats['total_runtime']:.2f}s")
        
        if 'rules' in stats:
            print(f"Rules executed: {len(stats['rules'])}")
            
            rule_times = []
            for rule, data in stats['rules'].items():
                if 'runtime' in data:
                    rule_times.append((rule, data['runtime']))
            
            if rule_times:
                rule_times.sort(key=lambda x: x[1], reverse=True)
                print("\nTop 5 slowest rules from stats:")
                for rule, runtime in rule_times[:5]:
                    print(f"  {rule}: {runtime:.2f}s")
    
    except Exception as e:
        print(f"Error reading stats file: {e}")

if __name__ == "__main__":
    import sys
    
    # Parse command line arguments
    benchmark_dir = "benchmarks"
    experiment = None
    
    if len(sys.argv) > 1:
        benchmark_dir = sys.argv[1]
    
    if len(sys.argv) > 2:
        experiment = sys.argv[2]
    
    # Try to auto-detect experiment from config if not provided
    if not experiment:
        try:
            import yaml
            with open("config.yaml", "r") as f:
                config = yaml.safe_load(f)
                experiment = config.get("experiment")
                if experiment:
                    print(f"🔍 Auto-detected experiment {experiment} from config.yaml")
        except:
            pass
    
    print("Starting benchmark analysis...")
    analyze_benchmarks(benchmark_dir, experiment)
    analyze_stats_file()
    print("\nAnalysis complete!")
