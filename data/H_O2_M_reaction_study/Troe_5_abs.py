Absorption-coefficients:
  - species: H2O2
    wave-lengths: 
      - value: 215
        unit: nm
        functional-form: A
        source: Hong, but fit to equation from Kappel
        parameter-one:
          value: 96000
          units: 'mol/c,'
          absolute-uncertainty: 
            value: .7
            units: 'cm^2/mol^-1'
        parameter-two:
          value: -19.0
          units: 'cm^2/K*mol^-1'
          absolute-uncertainty:
            value: .3
            units: 'cm^2/mol^-1'
  - species: HO2
    wave-lengths: 
      - value: 215
        unit: nm
        functional-form: B
        source: Hong, but fit to equation from Kappel
        parameter-one:
          value: 1155000
          units: 'cm^2/mol^-1'
          absolute-uncertainty:
            value: .7
            units: 'cm^2/mol^-1'
        parameter-two:
          value: -1299000
          units: 'cm^2/K*mol^-1'
          absolute-uncertainty:
            value: .3
            units: 'cm^2/K*mol^-1'