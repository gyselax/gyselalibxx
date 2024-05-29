# Reaction rates

The `reactionrates` folder contains all the code describing reactions within the `geometryMX` framework. 

The interface IReactionRates describes reaction rate that only depend on the temperature. Existing reaction rates that implements
such interface are : 
- ConstantIonizationRate
- ConstantRecombinationRate
- ConstantChargeExchangeRate