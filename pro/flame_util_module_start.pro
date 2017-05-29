PRO flame_util_module_start, fuel, module_name


	fuel.util.last_routine_time = systime(/seconds)
  fuel.util.last_routine_name = module_name


  print, ''
  print, '-------------------------------------'
  print, '---      ' + module_name + '        ---'
  print, '-------------------------------------'
  print, ''



END
