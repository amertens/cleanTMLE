# the printed design report is stable (snapshot)

    Code
      print(rep)
    Output
      Design report (pre-outcome)
        Roles:    analyst: programmer A; review_team: panel B
        Support: PASS; 0.0% outside [0.05, 0.95]; max weight 5; ESS 208 / 158
        Feasible estimands: ATE, trimmed_ATE, ATT, ATC, ATO, matched_ATT
        Declared ladder: ATE -> trimmed_ATE -> ATT (trigger SEVERE)
      
        The declared primary (ATE) is supported (overlap verdict PASS). Proceed to estimation.

# the Muntner-format decision log export is stable (snapshot)

    Code
      print(log, row.names = FALSE)
    Output
                      date                   stage
       <timestamp>         estimand_ladder
       <timestamp> Stage 1a (declarations)
       <timestamp>     Stage 4 (unmasking)
                                                                              issue
                         Primary ATE; fallbacks trimmed_ATE -> ATT; trigger SEVERE.
       Declared negative control(s): nc_outcome (domains: health_seeking_behavior).
                                        Outcome 'event_24' unmasked for estimation.
                              decision rationale  decided_by
                                  <NA>      <NA>        <NA>
                                  <NA>      <NA>        <NA>
       unmask and authorise estimation      <NA> review team

