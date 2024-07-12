# Reaction rates

The `reactionrates` folder contains all the code describing reactions within the `geometryMX` framework.

The interface IReactionRates describes reaction rate that depend on density and temperature.
Existing reaction rates that implements such interface are: 
- ConstantRate
- ChargeExchange
- Ionization
- Recombination