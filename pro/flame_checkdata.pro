PRO flame_checkdata, fuel

	flame_util_module_start, fuel, 'flame_checkdata'




  flame_util_module_end, fuel


  print, '-------------------------------------'
	print, '-------------------------------------'
	print, '-------------------------------------'

	; print total execution time
	print, ' '
	print, 'The data reduction took a total of ', $
		cgnumber_formatter((systime(/seconds) - fuel.util.start_time)/60.0, decimals=2), ' minutes.'
		print, ' '


END
