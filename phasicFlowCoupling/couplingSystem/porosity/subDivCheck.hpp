{
	realx3 subPoint = pPos + offset;

	// cell id of subPoint
	auto subCellId = cMesh_.findCellTree(subPoint, cellId);


	if(subCellId >= 0 )
	{
		if(cellId == subCellId )
		{
			numInCenter++;
		}
		else
		{
			solidVol[subCellId] += pSubVol;	
		}
		
	}

}