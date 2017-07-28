PRO flame_util_module_end, fuel

    filename = fuel.util.intermediate_dir + 'fuel.sav'
    save, fuel, filename=filename

    print, ''
    print, 'fuel structure saved to ' + filename
    print, fuel.util.last_routine_name + ' took ', $
      cgnumber_formatter( systime(/seconds) - fuel.util.last_routine_time, decimals=2), ' seconds'
    print, ''



END
