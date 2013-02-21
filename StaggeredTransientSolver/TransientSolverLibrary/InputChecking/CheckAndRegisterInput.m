function [] = CheckAndRegisterInput(~,Data)
     
%     ThrowingID = 'TransientSolverCheckAndRegisterInput:';
%     
% % ======================================================================== %
% %    Get the handles that generate the required/optional field arrays.     %
% % ======================================================================== %
%     switch(Switch)
%         case('SolverSetup'),    GetFields = @FieldsForSolverSetup;
%         case('System'),         GetFields = @FieldsForSystem;
%         case('ControlVolumes'), GetFields = @FieldsForControlVolumes;
%         case('Interfaces'),     GetFields = @FieldsForInterfaces;
%         case('Miscellaneous'),  GetFields = @FieldsForMiscellaneous;
%             
%         otherwise
%             error([ThrowingID,'ImproperSwitch'],...
%                 'The switch case ''%s'' is not valid.',Switch);
%     end
%     [Required,Optional,Default,Shorten,ShortenMask] = GetFields();
    
    
% ======================================================================== %
%                Get the field names and proccess the data                 %
% ======================================================================== %
    FieldNames = fieldnames(Data);
    
    for k = 1:length(FieldNames)
        
        FieldName = FieldNames{k};
        FieldData = Data.(FieldName);
        Database('insert',FieldName,FieldData);

    end
    
    
end

% function [Required,Optional,Default,Shorten,ShortenMask] = SolverSetupFields()
%     Required    = {'TimeStart','TimeEnd' };
%     Optional    = {'dTimeSave','CFLLimit'};
%     Default     = {     1     ,   0.8    };
%     Shorten     = {};
%     ShortenMask = {};
% end

