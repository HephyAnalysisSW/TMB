# create event-wise coeffs with photon_pt, adapted from WeightInfo.py
def get_coeff_list_with_photon_pt_from_events(sample, selectionString=None, weightFunction=None):
        ''' Create list of weights for each event
        '''
        # RootTools
        from RootTools.core.standard             import TreeVariable, VectorTreeVariable

        sample.setSelectionString(selectionString) 
        
        variables = map( TreeVariable.fromString, [ "np/I", "photon_pt/F" ] )
        variables.append( VectorTreeVariable.fromString('p[C/F]', nMax=1000) )

        reader = sample.treeReader( variables = variables )
        reader.start()

        coeffs_with_photon_pt = []
        while reader.run():
            event_data = {
                'coeffs': [reader.event.p_C[i] * (weightFunction( reader.event, sample )) if weightFunction is not None else reader.event.p_C[i] for i in range(reader.event.np)],
                'photon_pt': reader.event.photon_pt
            }
            coeffs_with_photon_pt.append(event_data)

        return coeffs_with_photon_pt