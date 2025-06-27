  using Pkg
  Pkg.add("Coverage")

  using Coverage

  # Process coverage data
  coverage = process_folder()

  # Get coverage summary
  covered_lines, total_lines = get_summary(coverage)
  println("Coverage: $(round(100 * covered_lines / total_lines, digits=2))%")

  # Generate LCOV format (what CodeCov expects)
  LCOV.writefile("lcov.info", coverage)

    # Install lcov tools if you haven't already
  brew install lcov  # on macOS

  # Generate HTML report
  genhtml lcov.info -o coverage_report

  # Open in browser
  open coverage_report/index.html



  # Run tests with coverage
  run(`julia --project=. --code-coverage=user -e 'using Pkg; Pkg.test()'`)

  # Process and generate LCOV
  coverage = process_folder()
  LCOV.writefile("lcov.info", coverage)

  # Show summary
  covered, total = get_summary(coverage)
  println("\nCoverage: $(round(100covered/total, digits=2))%")

  # Clean up
  clean_folder(".")

  # Generate HTML
  run(`genhtml lcov.info -o coverage_report`)
  println("\nHTML report generated in coverage_report/index.html")